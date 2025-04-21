clc
clear all
close all
%%
load('PFCspikeTimes.mat') % spikeTimes in ms
spikeIntervals = spikeTimes(2:length(spikeTimes)) - spikeTimes(1:length(spikeTimes) - 1);% inter-spike interval
figure(1);                                              % use Figure 1
binSize = 1;                                            % 1 ms bins
x = 1:binSize:200;
intervalDist = hist(spikeIntervals(spikeIntervals < 200), x);
intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number
bar(x, intervalDist);
set(gca,'FontSize',20,'XLim',[0 200]); xlabel('ISI (ms)'); ylabel('count');
%%
% As promised, the histogram of interspike intervals looks like an exponential 
% probability distribution.  For verificiation, we'll fit it with an exponential
% probability density function.
[maxISI,maxid]= max(intervalDist);
xdata = maxid:binSize:200;
intervalDist_decay = intervalDist(maxid:end);
start_point = 20;
expfun = @(b,xdata)  1/b * exp(-xdata/b);  
exponentialfit = @(b) sum((intervalDist_decay - expfun(b,xdata)).^2);   
[lambda_expest, SSE_exp] = fminsearch(exponentialfit, start_point);
y_expfit = expfun(lambda_expest,xdata);            % exponential function
hold on;
plot(xdata, y_expfit, 'r:','linewidth',2);      
legend1 = ['Exponential, ', 'SSE=',num2str(SSE_exp)];
%%
% Let's compare it by fitting with a linear
% probability density function.
start_point = [0.2,0.02];
linear_fun = @(b,xdata)  (-b(1) * xdata + b(2));  
linearfit = @(b) sum((intervalDist_decay - linear_fun(b,xdata)).^2);   
[lambda_linest, SSE_lin] = fminsearch(linearfit, start_point);
y_linfit = linear_fun(lambda_linest,xdata);            % linear function
xlabel('Interspike interval');
ylabel('Probability');
plot(xdata, y_linfit, 'y:','linewidth',2);  
legend2 = ['Linear, ', 'SSE=',num2str(SSE_lin)];
%%
% Let's compare it by fitting with an inverse
% probability density function.
start_point = 1;
inverse_fun = @(b,xdata) b(1)./(xdata);
inversefit = @(b) sum((intervalDist_decay - inverse_fun(b,xdata)).^2);   
[lambda_invest, SSE_inv] = fminsearch(inversefit, start_point);
y_invfit = inverse_fun(lambda_invest,xdata);            % inverse function
xlabel('Interspike interval');
ylabel('Probability');
hold on;
plot(xdata, y_invfit, 'g:','linewidth',2);                                      
hold off;
axis([min(x) max(x) 0 maxISI * 1.1]); 
legend3 = ['Inverse, ', 'SSE=',num2str(SSE_inv)];
legend(legend1,legend2,legend3);
