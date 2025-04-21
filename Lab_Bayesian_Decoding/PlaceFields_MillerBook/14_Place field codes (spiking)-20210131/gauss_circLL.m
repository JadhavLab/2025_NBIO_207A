function nLL=gauss_circLL(params, positions, spike_vector)

% Takes an array of position coordinates combined with the spike train
% synchronized in time and calculates the log likelihood of the spikes 
% being produced by a given Gaussian receptive field. 
%
% In this model the amplitude parameter will be equivalent to rmax*dt,
% which is the maximum probability of a spike at the center of the
% Gaussian.

Inputx = positions(:,1);       % x locations of spikes
Inputy = positions(:,2);       % y locations of spikes

A=params(1);                    % spike rate*dt at center of place field
x0=params(2);                   % x-coordinate of place field center
y0=params(3);                   % y-coordinate of place field center
sigma=params(4);                % standard deviation of place field
xdev = Inputx - x0;             % x-positions relative to place-field center
ydev = Inputy - y0;             % y-positions relative to place-field center

spike_times = find(spike_vector);   % time bins when the cell emits a spike

% xs and ys are distances in x- and y-directions respectively from the 
% center of the receptive field of positions of spikes
xs = xdev(spike_times);
ys = ydev(spike_times);

% log likelihood of all the spikes at the positions they occurred
ll_spikes = sum(log(A)-  (xs.*xs + ys.*ys) /(2*sigma*sigma)  );


nspike_times = find(~spike_vector); % time bins without the cell spiking

% x- and y- distances from center of reeptive field at times of no spikes
xns = xdev(nspike_times);
yns = ydev(nspike_times);

% log likelihood sum over all time bins without spikes
ll_nspikes = sum( log( 1-A*exp(-(xns.*xns +yns.*yns)/(2*sigma*sigma)) ) );

% total log-likelihood is sum of contributions from times with spikes and
% times without spikes
LL = ll_spikes + ll_nspikes ;

% Return the negative of the log likelihood as the function will be
% minimized
nLL = -LL;

% meanx and meany denote x,y coordinates of the center of the environment
meanx = (min(Inputx)+max(Inputx))/2;
meany = (min(Inputy)+max(Inputy))/2;    

% xrange and yrange are distances from the center to the edges of the
% environment in the x- and y- directions respectively
xrange = (-min(Inputx)+max(Inputx))/2;
yrange = (-min(Inputy)+max(Inputy))/2;

% xoffset and yoffset indicate the distances of the estimated receptive
% field center from the center of the environment in the x- and
% y-directions respectively
xoffset = x0-meanx;
yoffset = y0 - meany;

% In the following section penalize sigma if it is too large and similarly
% the xoffset and yoffset if they are too far from the center of the grid
% (i.e. outside the range of the environment).
nLL = nLL + 0.1*length(Inputx)*(sigma^2/(xrange^2+yrange^2+sigma^2) ) + ...
    0.1*length(Inputx)*(xoffset*xoffset/(xrange*xrange+xoffset*xoffset) )+ ...
    0.1*length(Inputx)*(yoffset*yoffset/(yrange*yrange+yoffset*yoffset) );