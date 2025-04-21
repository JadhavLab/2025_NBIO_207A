clc
clear all
close all;
%%
% load data
load('PlaceCells_demo_data.mat')
tBinSz_sm = 200; % decoding window tau, 200ms
fx_placefield = fx_placefield+ (eps.^8); %f(x), Add a small number so there are no zeros
maxcellnum = length(fx_placefield(1,:)); % total number of cells, N
timevec = 0:tBinSz_sm/1000:max(truex_position(:,1));% time vector
%%
% let's start with just one cell
cellnum = 1; % number of cells, N
cellID = randi(maxcellnum,cellnum); %randomly pick one cell
fx_thiscell = fx_placefield(:,cellID); %tuning curve of this cell
y_thiscell = y_spiketrain(find(y_spiketrain(:,1) == cellID),:);% spike train of this cell
expecSpk = fx_thiscell.*tBinSz_sm./1000; %Expected number of spikes per bin,i.e. tau*f(x)

%%
% let's what the cell response looks like
figure,
subplot(211)
plot(x_position, fx_thiscell)
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')
title('Tuning curve (place field)')
subplot(212)
running_ID = find(y_thiscell(:,1)>0);
plot(y_thiscell(running_ID,2),y_thiscell(running_ID,1),'.')
xlim([0,max(timevec)])
ylim([cellID-0.5,cellID+0.5])
xlabel('Time (s)')
ylabel('Cell No.')
title('Spike train')
%%
% let's do bin-by-bin decoding using the spike train from the cell
pMat_all = [];% to save the posterior probablity matrix
currentpos_all = [];% to save the current positions

% time bin loop
for i = 1:length(timevec)
    bin = timevec(i);% current time
    [~,avglindist_bin] = min(abs(truex_position(:,1) - bin));% current position ID
    avglindist = truex_position(avglindist_bin(1),2);% current position
    currentpos_all = [currentpos_all,avglindist];% save current postion

    spikebins = y_thiscell(find(y_thiscell >= (bin - tBinSz_sm/2000) & y_thiscell <= (bin + tBinSz_sm/2000))); % spikes in the current bin
    spikecount = length(spikebins); % spike count
    if (spikecount > 0) % decode the bin that has spikes
        nPBin = size(fx_thiscell); % N position bin
        wrking = fx_thiscell .^ spikecount; %f(x)^y
        expon = exp(-expecSpk); %e^-(tau*f(x))
        pMat = wrking.*expon; %f(x)^y * e^-(tau*f(x))
        pMat(isnan(pMat)) = 0; % set empty bins to zero
        pMat = pMat./sum(pMat);% C(x), normalization to make sure sum(p(x|y))=1
        %hence, pMat is the posterior probability matrix p(x|y) = C*prod(f(x)^y) * e^-sum(tau*f(x))
    else
        pMat = zeros(size(fx_thiscell));
    end
    pMat_all = [pMat_all,pMat];% pMat is the decoding result for each bin, and pooled into pMat_all
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
running_ID = find(y_thiscell(:,1)>0);
plot(y_thiscell(running_ID,2),y_thiscell(running_ID,1),'.')
xlim([0,max(timevec)])
ylim([cellID-0.5,cellID+0.5])
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
