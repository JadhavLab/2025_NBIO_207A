cd('F:\Research\Code\2025_NBIO_207A\Lab_SpikeTrains_PoissonProcesses\Matlab_SpikeData_PointProcess_ProblemSet')
load('hipp_data.mat')

%% 1 Plot rat in 2D trajectory with spikes from Neuron 1 and 2 overlaid
spikes = logical(spikes);
spikes2 = logical(spikes2);

figure; hold on; axis square
plot(xN,yN,'Color',[.5 .5 .5])
plot(xN(spikes),yN(spikes),'r.')
plot(xN(spikes2),yN(spikes2),'b.')

%% 2 Fit a GLM with linear position to this point process data

predictors1 = [xN yN];
response1 = spikes;
[model1_b,model1_dev,model1_stats] = glmfit(predictors1,response1,'poisson');
% model1_b = vector of coefficient estimates by MLE (default)
% model1_dev = -2logL(theta) deviance of the model fit
% model1_stats = struct of statistics on the fitted model

model1_aic = model1_dev + 2*(size(predictors1,2)) + 1;

% Define custom link function 'poisson' for glmval (not builtin)
% link1 = @(mu) poissinv(mu);
% derlink1 = @(mu) 1 ./ poisspdf(poissinv(mu));
% invlink1 = @(resp) cdf('poisson',resp);
% F = {link1, derlink1, invlink1};
[spikes_pred1,ciHi1,ciLo1] = glmval(model1_b, predictors1,'logit',model1_stats);
% predY = exp(model1_b(1) + model1_b(2)*xN + model1_b(3)*yN);

figure; hold on
plot(spikes.*max(spikes_pred1),'k')
plot(spikes_pred1,'b'); 
plot(spikes_pred1+ciHi1,'Color',[.5 .5 .5]); plot(spikes_pred1-ciLo1,'Color',[.5 .5 .5])

%% 3 Fit GLM with quadratic position information

predictors2 = [xN xN.^2 yN yN.^2];
response2 = spikes;
[model2_b,model2_dev,model2_stats] = glmfit(predictors2,response2,'poisson');

model2_aic = model2_dev + 2*(size(predictors2,2)) + 1;

if model2_aic < model1_aic
    disp(['Model 2 AIC lower (better) than Model 1'])
end

[spikes_pred2,ciHi2,ciLo2] = glmval(model2_b, predictors2,'log',model2_stats);

figure; hold on
plot(spikes.*max(spikes_pred2),'k')
plot(spikes_pred1,'b'); 
plot(spikes_pred2,'r'); 
plot(spikes_pred1+ciHi1,'Color',[.5 .5 .5]); plot(spikes_pred1-ciLo1,'Color',[.5 .5 .5])
plot(spikes_pred2+ciHi2,'Color',[.5 .5 .5]); plot(spikes_pred2-ciLo2,'Color',[.5 .5 .5])
% plot(smooth(spikes,250))

%% Create a simulated spike train from the predicted lambdas - stochastic
% May be useful for visualization purposes or comparing true spike train to
% inferred spike train of GLM model

spike_trn2 = poissrnd(spikes_pred2);    %Treat predicted spike rate as lambda

figure; hold on;
plot(spikes)
plot(spike_trn2)







