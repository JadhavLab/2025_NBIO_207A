clc
clear all
close all;
%%
% load data
load('PlaceCells_demo_data.mat')
tBinSz_sm = 200; % decoding window tau, 200ms
fx_placefield = fx_placefield+ (eps.^8); %f(x), Add a small number so there are no zeros
maxcellnum = length(fx_placefield(1,:)); % number of cellsm, N
timevec = 0:tBinSz_sm/1000:max(truex_position(:,1));% time vector
%%
% let's try 5 cells this time
cellnum = 5; % number of cells, N
cellID = randi(maxcellnum,cellnum); %randomly pick 5 cells
fx_placefield_5cells = fx_placefield(:,cellID); %tuning curve of these 5 cells
expecSpk = fx_placefield_5cells.*tBinSz_sm./1000; %Expected number of spikes per bin,i.e. tau*f(x)

% we sort these 5 cells first
[cellID_new,index] = sort(cellID);
fx_placefield_5cells = fx_placefield_5cells(:,index);
y_spiketrain_5cells = [];
for i = 1:cellnum
    thiscell_spikes = y_spiketrain(find(y_spiketrain(:,1) == cellID(index(i))),2);% spike train of this cell
    y_spiketrain_5cells = [y_spiketrain_5cells;i*ones(size(thiscell_spikes)),thiscell_spikes];
end
%%
% let's what the response of these 5 cells look like
figure,
subplot(211)
plot(x_position, fx_placefield_5cells)
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')
title('Tuning curve (place field)')
subplot(212)
running_ID = find(y_spiketrain_5cells(:,1)>0);
plot(y_spiketrain_5cells(running_ID,2),y_spiketrain_5cells(running_ID,1),'.')
xlim([0,max(timevec)])
ylim([0,5])
xlabel('Time (s)')
ylabel('Cell No.')
title('Spike train')
%%
% let's do bin-by-bin decoding using the spike train from the 5 cells

pMat_all = [];% to save the posterior probablity matrix
currentpos_all = [];% to save the current positions

% time bin loop
for i = 1:length(timevec)
    bin = timevec(i);% current time
    [~,avglindist_bin] = min(abs(truex_position(:,1) - bin));% current position ID
    avglindist = truex_position(avglindist_bin(1),2);% current position
   
    % get all the spikes in the current bin
    clear spikecount
    celldata = [];
    for n = 1:cellnum
        spikes_i  = y_spiketrain_5cells(find(y_spiketrain_5cells(:,1) == n),2);%spikes of the current cell
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
        
        fx_placefield_active = fx_placefield_5cells(:,cellsi);% f(x)
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
    currentpos_all = [currentpos_all,avglindist];% save current postion
end
%%
% visualization of the decoding result
figure,
subplot(211)
imagesc(-1.*pMat_all) %reverse the axis
ax = gca;
ax.YDir = 'normal';
hold on
colormap(gray)
plot(currentpos_all/2,':','color','cyan','linewidth',2)
caxis([-0.3,0])%probability 0-0.3
axis off
title('Decoding result')
subplot(212)
running_ID = find(y_spiketrain_5cells(:,1)>0);
plot(y_spiketrain_5cells(running_ID,2),y_spiketrain_5cells(running_ID,1),'.')
xlim([0,max(timevec)])
ylim([0,5])
xlabel('Time (s)')
ylabel('Cell No.')
title('Spike train')
%%
% visualization of the decoding error
[~,estimated_posID] = max(pMat_all);
estimated_pos = x_position(estimated_posID)';
decode_error = abs(estimated_pos-currentpos_all);
figure,
counts = hist(decode_error,0:0.5:200);
plot(0:0.5:200,cumsum(counts)./sum(counts),'k','linewidth',2) % cumulative histogram
hold on
plot(median(decode_error),0.5,'r.','markersize',40) % median decoding error
hold off
xlabel('Decoding error (cm)')
ylabel('Cumulative fraction')
title(['Decoding error, median = ',num2str(round(median(decode_error)*100)/100),'cm'])
