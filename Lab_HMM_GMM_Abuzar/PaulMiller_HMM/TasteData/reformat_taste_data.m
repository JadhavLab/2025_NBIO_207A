% reformat_taste_data.m
%
%   This extracts responses after the 2000ms time point of stimulus onset
%   and can produce multiple time bins, which produces multiple responses
%   per neuron per taste in a manner that increases the dimensionality of
%   the array holding the data.
%
%   Also, the code removes trials in which the number of spikes is well
%   below the mean, as there are known to be faulty trials in the data set.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [analysis_data, id_data] = reformat_taste_data(ps_array)

% Extract from the array the numbers of stimuli, trials, cells, and time
% points
[ntastes, ntrials, ncells, nt] = size(ps_array);    

% Following lines determine how to split the activity into time bins
% To just use total spike counts, set Nbins to 1 and binwidth to a large
% number of milliseconds like 2500.
Nbins = 1;             % Number of distinct time bins
binwidth = 2500;         % Width of each time bin
binsep = 100;           % Separation of onset times of time bins

% counts is an array where each trial is a separate row and each column is
% the number of spikes of a particular cell in a particular time bin
counts = zeros(ntrials*ntastes,ncells*Nbins);   
taste_ids = zeros(ntrials*ntastes,1);       % Labels the taste id of each row

% Now loop through time bins and fill up the array
for i = 1:Nbins
    binstart = 2001 + (i-1)*binsep;         % bin onset (stimulus onset is 2000ms)
    binstop = binstart + binwidth-1;        % end of time bin
    
    % In the next line count spikes (sum) in all 1ms-bins across the time
    % bin and reduce the rank of the array
    total_ps_array = squeeze(sum(ps_array(:,:,:,binstart:binstop),4));
    
    % The data goes into a set of columns that are cell specific but
    % shifted by the number of time bins multiplied by number of cells
    cellstart = (i-1)*ncells + 1;
    cellstop = i*ncells;

    % Now put data for all trials of each taste stimulus in a different set
    % of rows from the "total_ps_array" (with 3 indices) into the counts
    % array (which has 2 indices), again using squeeze.
    for taste = 1:ntastes
        counts(ntrials*(taste-1)+1:ntrials*taste,cellstart:cellstop) = ...
            squeeze(total_ps_array(taste,:,:));
        taste_ids(ntrials*(taste-1)+1:ntrials*taste) = taste;
        
    end
end

% The next section finds the number of spikes in each trial and removes
% trials with too few spikes.
sum_counts = sum(counts,2);             % No. of spikes in each trial
mean_sum = mean(sum_counts);            % Mean no. of spikes
bad_indices = find(sum_counts < mean_sum/10);   % Trials with < 10% of mean
num_bad = length(bad_indices);                  % No. of such trials
disp(strcat('numbad = ',num2str(num_bad)));     % Display no. of trials we will remove

trials_to_include = ones(1,ntrials*ntastes);    % By default include all trials
trials_to_include(bad_indices) = 0;             % Set to 0 for trials with few spikes
includes = find(trials_to_include);             % Only include the trials with 1
analysis_data = counts(includes,:);             % Array with included trials only
id_data = taste_ids(includes);                  % Array with taste ids of included trials

