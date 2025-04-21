% Fourier_gaussian.m
%
% Combines a Gaussian with the Fourier transform then shifts through time
% to produce the spectral density
clear

tic
% tvec = 0.0001:0.0001:1;
% data = cos(2*pi*40*tvec);
% data(round(end/2):end) = cos(2*pi*40*tvec(round(end/2):end)+pi/2);

load('V_coarse.mat')
dt = 0.001;
data = V_coarse(1:end/1000) - mean(V_coarse(1:end/1000));
clear V_coarse;

Nt = length(data);
tmax = dt*Nt;

Gauss_sig = 5;
n_sig = round(Gauss_sig/dt);
two_sig_sq = 2*n_sig*n_sig;

f_vec = 0:0.2:300;          % Set of frequencies to test
Nf = length(f_vec);         % Number of distinct frequencies

i_tvec = [1:Nt]';
coarse_i_tvec = 50:50:Nt;
tvec = i_tvec*dt;
index = 0;
for i_t = coarse_i_tvec;
    i_t
    index = index + 1;
    gauss_curve = exp(-(i_tvec-i_t).*(i_tvec-i_t)/two_sig_sq);  % Gaussian window
    
    new_data = data.*gauss_curve;       % Data multiplied by the Gaussian window
    
    for i_f = 1:Nf;               % Loop through set of frequencies
        f = f_vec(i_f);           % Value of frequency to test
        sinvec = sin(2*pi*f*tvec);      % sine function at that frequency
        cosvec = cos(2*pi*f*tvec);      % cosine function at that frequency
        
        % Now take the normalized dot produce of the sine function then the
        % cosine function with the firing rate vector.
        sin_overlap = mean(sinvec.*new_data);
        cos_overlap = mean(cosvec.*new_data);
        
        A(index,i_f) = sin_overlap;
        B(index,i_f) = cos_overlap;
    end
    
end

power = A.*A + B.*B;


imagesc(power')

set(gca,'YDir','Normal')

theta = atan(A./B) +pi/2 * (sign(B)-1);
figure()
imagesc(theta')
set(gca,'YDir','Normal')
toc
