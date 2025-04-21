function nLL=gauss2DLL(params, positions, spike_vector)

% Takes an array of position coordinates combined with the spike train
% synchronized in time and calculates the log likelihood of the spikes 
% being produced by a given Gaussian receptive field. 
%
% In this model the amplitude parameter will be equivalent to rmax*dt,
% which is the maximum probability of a spike at the center of the
% Gaussian.

Inputx = positions(:,1);       % x locations of spikes
Inputy = positions(:,2);       % y locations of spikes

A=params(1);                    % Amplitude, r*dt at peak rate, r
x0=params(2);                   % x-coordinate of center
y0=params(3);                   % y-coordinate of center
sigma1=params(4);               % s.d. along 1st axies
sigma2 = params(5);             % s.d. along 2nd axis
theta = params(6);              % angle of rotation between x- and 1st-axis

xdev = Inputx - x0;             % x-axis distance from spikes to place field center
ydev = Inputy - y0;             % y-axis distance from spikes to place field center
x1 = xdev*cos(theta) + ydev*sin(theta);     % 1st axis distance from spikes to place field center
x2 = -xdev*sin(theta) + ydev*cos(theta);    % 2nd axis distance from spikes to place field center

spike_times = find(spike_vector);           % Time bins containing spikes

x1s = x1(spike_times);          % 1st axis distance from spike positions to place field center
x2s = x2(spike_times);          % 2nd axis distance from spike positions to place field center

% sum of log-likelihood of all spike probabilities given the place field.
ll_spikes = sum(log(A)-(x1s.*x1s)/(2*sigma1*sigma1) -(x2s.*x2s)/(2*sigma2*sigma2) );

nspike_times = find(~spike_vector);      % Time bins containing no spikes

x1ns = x1(nspike_times);    % 1st axis distance from "no-spike" positions to place field center
x2ns = x2(nspike_times);    % 2nd axis distance from "no-spike" positions to place field center

% sum of log-likelihood of no spike in all time bins without a spike 
ll_nspikes = sum(log(1-A*exp(-(x1ns.*x1ns)/(2*sigma1*sigma1) -(x2ns.*x2ns)/(2*sigma2*sigma2)) ) );

% Total log likelihood is some of time bins with spikes plus time bins
% without spikes
LL = ll_spikes + ll_nspikes;

nLL = -LL;  % Return negative of log-likelihood, which must be minimized
