clc
clear all
close all;
%%
% load data
load('PlaceCells_demo_data.mat')
tBinSz_sm = 20; % decoding window tau, 20ms
fx_placefield = fx_placefield+ (eps.^8); %f(x), Add a small number so there are no zeros
cellnum = length(fx_placefield(1,:)); % number of cellsm, N
expecSpk = fx_placefield.*tBinSz_sm./1000; %Expected number of spikes per bin,i.e. tau*f(x)
load('ReplayEvent2.mat')
%%
% bin-by-bin decoding
timevec = replay_time(1):tBinSz_sm/1000:replay_time(2);% time vector
pMat_all = [];% to save the posterior probablity matrix

% time bin loop
for i = 1:length(timevec)
    bin = timevec(i);% current time
    
    % get all the spikes in the current bin
    clear spikecount
    celldata = [];
    for n = 1:cellnum
        spikes_i  = y_spiketrain(find(y_spiketrain(:,1) == n),2);%spikes of the current cell
        if ~isempty(spikes_i)
            spikebins = spikes_i(find(spikes_i >= (bin - tBinSz_sm/2000) & spikes_i <= (bin + tBinSz_sm/2000)));
            spikecount(n) = length(spikebins);
            tmpcelldata = [length(spikebins),n];
        else
            tmpcelldata = [0,n];
            spikecount(n) = 0;
        end
        celldata = [celldata;tmpcelldata];
    end
    
    cellcounts = sum((spikecount > 0));
    if (cellcounts > 0) % decode the bin that has spikes
        cellsi = find(spikecount > 0);%active cell
        spkPerBin = celldata(cellsi,1)'; % spike counts per bin
        
        fx_placefield_active = fx_placefield(:,cellsi);% f(x)
        nPBin = size(fx_placefield_active,1); % N position bin
        wrking = bsxfun(@power, fx_placefield_active, spkPerBin); %f(x)^y
        wrking = prod(wrking,2);%prod(f(x)^y)
        expon = exp(-sum(expecSpk,2)); %e^-sum(tau*f(x))
        pMat = bsxfun(@times,wrking, expon); %prod(f(x)^y) * e^-sum(tau*f(x))
        pMat(isnan(pMat)) = 0; % set empty bins to zero
        pMat = pMat./sum(pMat);% C(x), normalization to make sure sum(p(x|y))=1
        %hence, pMat is the posterior probability matrix p(x|y) = C*prod(f(x)^y) * e^-sum(tau*f(x))
    else
        pMat = zeros(size(x_position));
    end
    pMat_all = [pMat_all,pMat];% pMat is the decoding result for each bin, and pooled into pMat_all
end
%%
% visualization of the decoding result
figure,
subplot(211)
imagesc(pMat_all) %reverse the axis
ax = gca;
ax.YDir = 'normal';
hold on
colormap(hot)
caxis([0,0.1])%probability 0-0.1
axis off
subplot(212)
replay_ID = find(y_spiketrain(:,2)>= replay_time(1) & y_spiketrain(:,2)<= replay_time(2));
plot(y_spiketrain(replay_ID,2)-replay_time(1),y_spiketrain(replay_ID,1),'.')
xlim([0,replay_time(2)-replay_time(1)])
ylim([0,cellnum])
xlabel('Time (s)')
ylabel('Cell No.')
