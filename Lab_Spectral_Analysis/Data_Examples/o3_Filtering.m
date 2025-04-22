%% Filtering - Signal Processing Toolbox
% -----------------------------------------

% Difference equation: What is it?

%example - running avaerage
% http://www.mathworks.com/help/matlab/data_analysis/filtering-data.html

load count.dat
x = count(:,1);
t = 1:length(x);
a = 1; % a_0 is the (hidden) coefficient on the left side, in front of y(n)
b = [1/4 1/4 1/4 1/4]; % four b's of 1/4 each so we get the mean
y = filter(b,a,x); % x is the original signal, y the filtered version

figure; hold on;
plot(t,x,'-.',t,y,'-'), grid on
legend('Original','Filtered',2)


keyboard;


%% Filters- high pass/ low pass/ band pass
% eg. Butterworth vs. Chebyshev

help butter
% butter Butterworth digital and analog filter design.
% [B,A] = butter(N,Wn) designs an Nth order lowpass digital
% Butterworth filter and returns the filter coefficients in length
% N+1 vectors B (numerator) and A (denominator). The coefficients
% are listed in descending powers of z. The cutoff frequency
% Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding to
% half the sample rate.
% If Wn is a two-element vector, Wn = [W1 W2], butter returns an
% order 2N bandpass filter with passband W1 < W < W2.


%Generate a 10-second long white noise signal, sampled at 500Hz. Filter it using a bandpass Butterworth filter between 50 and 100 Hz of order 4.
%Plot the Welch spectrum, in dB, of the original and filtered signal, using a 512-sample Hanning window. Evaluate the FFT over 2^14 points.

% set up time axis
Fs = ...
tvec = ...
% generate white noise
x = rand(...)
% get PSD
[Porig,Forig] = pwelch(x, ...)
% design filter
W1 = ...
W2 = ...
[b,a] = butter(4,[W1 W2]);
y = filter(...)
% get PSD
[Pfilt,Ffilt] = pwelch(y, ...)
% plot the resulting PSDs
figire; hold on;
subplot(121)
plot(... 10*log10(..));
subplot(122)


% Determine filter order to pass as argument
Wp = [ 50 100] * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ 45 105] * 2 / Fs; % stopband
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[b2,a2] = butter(N,Wn); % builds filter
fvtool(b,a,b2,a2)

% Chebyshev

Wp = [ 50 100] * 2 / Fs;
Ws = [ 48 102] * 2 / Fs;
[N,Wn] = cheb1ord( Wp, Ws, 3, 20);
[b_c1,a_c1] = cheby1(N,0.5,Wn);
fvtool(b2,a2,b_c1,a_c1)

keyboard;

%% Phase responses and filtfilt
% ------------------------------

Fs = 500; dt = 1./Fs;
t = [0 10];
tvec = t(1):dt:t(2)-dt;
s1 = sin(2*pi*80*tvec+pi/6);
s2 = sin(2*pi*40*tvec);
s = s1 + s2;
sf = filter(b_c1,a_c1,s);

plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);

keybaord;

% Phase shifts
% Run fvtool again on the Butterworth and Chebyshev filters above, 
% and now select the Phase Response button in the top left of the window.

% Run filter forward and backward to eliminate phase shifts
sf = filtfilt(b_c1,a_c1,s);
plot(tvec,s,'k',tvec,sf,'r--'); hold on;
legend({'original','filtered'});
xlim([0 0.2]);


%% compare freq responses
Fs = 500; dt = 1./Fs;

t = [0 10];
tvec = t(1):dt:t(2)-dt;
x = rand(size(tvec)); % white noise input
[P,F] = pwelch(x,hanning(512),256,2^14,Fs);
y1 = filter(b_c1,a_c1,x);
[P1,F1] = pwelch(y1,hanning(512),256,2^14,Fs);
y2 = filtfilt(b_c1,a_c1,x);
[P2,F2] = pwelch(y2,hanning(512),256,2^14,Fs);
plot(F,10*log10(P),F,10*log10(P1),F,10*log10(
legend({'original','filter','filtfilt'});

keyboard;


%% Notch filter to remove noise

[b,a] = butter(10, [59 61] * 2 / Fs, 'stop');
fvtool(b,a);

% Chebyshev, or second-order section format
[z,p,k] = butter(10, [59 61] * 2 / Fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object
fvtool(h);


