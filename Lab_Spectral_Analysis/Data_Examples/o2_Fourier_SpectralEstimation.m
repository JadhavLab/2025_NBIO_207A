% Fourier series, transforms, power spectra
% -----------------------------------------

% See also MAtlab documentation on Spectral Analysis
% 

set(0,'DefaultAxesFontSize',18)


%% Sinusoids - Periodic functions of time / oscillations
% plot a simple sinusoid
Fs = 100; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1./Fs:t1; % construct time axis

f = 2; % frequency of sine to plot
y = sin(2*pi*f*tvec); % note sin() expects arguments in radians, not degrees (see also ''sind()'')
figure; hold on;
stem(tvec,y);

keyboard;
% Phase and amplitude of sinusoids
% -------------------------------
figure; hold on;
phi = pi/2;
subplot(221)
y = sin(2*pi*f*tvec + phi); % a phase shift
stem(tvec,y);
hold on;
plot(tvec,cos(2*pi*f*tvec),'r--','LineWidth',2) % notice, cosine is simply phase shifted sine
legend('sin (phase-shifted)', 'cos');
a = 2;
subplot(222)
y = a.*sin(2*pi*f*tvec + phi); % amplitude change
stem(tvec,y); % note scale of y axis!

keyboard;


% Frequency modulation of carrier signal: similar to cross-frequency coupling
% ------------------------------------------------------
figure; hold on;
f2 = 10;
m = 2;
subplot(311)
s1 = sin(2*pi*f*tvec);
plot(tvec,s1); title('message');
subplot(312);
s2 = sin(2*pi*f2*tvec);
plot(tvec,s2); title('carrier');
subplot(313);
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2));
plot(tvec,s3); title('FM signal');

keyboard;


%% Fourier analysis. Sum of sinusoids and harmonic series
% ------------------------------------------------------

%% harmonic series example
figure; hold on;
mag = [0.1 0 1.3 0.5]; % magnitudes for each term
pha = [-pi/6 0 pi 2*pi/3]; % phases for each term
f = 2; % base frequency
signal_out = zeros(size(tvec));
for ii = 1:numel(mag) % dont use i for index!
    this_signal = mag(ii)*cos(2*pi*f*ii*tvec + pha(ii));
    plot(tvec,this_signal,'r:'); hold on;
    signal_out = signal_out + this_signal; % build the sum
end
plot(tvec,signal_out,'LineWidth',2);
title('Signal created by sinusoid series')


keyboard;

% Fourier: any continuous signal can be represented by a series of sinusoids



%% Decomposing and reconstructing a signal
% ------------------------------------------------------
% FFT  computes the DFT to get amplitudes and phases of a signal
%%
rng('default'); % reset random number generator to reproducible state, so your plot will look like mine!
x = round(rand(1,8)*10); % generate a length 8 vector of integers between 0 and 10
xlen = length(x);
% get magnitudes and phases of Fourier series
X = fft(x);
Xmag = abs(X); % magnitudes, a_n
Xphase = angle(X); % phases, phi_n
n = 0:xlen-1;
t = 0:0.05:xlen-1; % a finer timescale to show the smooth signal later
for iH = xlen-1:-1:0 % reconstruct each harmonic
    s(iH+1,:) = Xmag(iH+1)*cos(2*pi*n*iH/xlen + Xphase(iH+1))/xlen;
    sm(iH+1,:) = Xmag(iH+1)*cos(2*pi*t*iH/xlen + Xphase(iH+1))/xlen;
    % detail: xlen appears here because the fundamental frequency used by fft() depends on this
end
ssum = sum(s); % coarse timescale (original points)
smsum = sum(sm); % fine timescale (to see full signal)
figure;
plot(n, x, 'go', t, smsum, 'b', n, ssum, 'r*');
legend({'original','sum - all','sum - points only'});

keyboard;


%% FFT output - magnitude and phase spectra

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1/Fs:t1-(1/Fs); % construct time axis; generate exactly 20 samples
f = 2; % signal frequency
y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s
yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
figure;
stem(yfft_mag);

% What is the frequency axis?   0-Fs/2
% Both real and complex components are returned

% Single-sided spectrum
Npoints = length(y);
F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis
yfft_mag = fftshift(yfft_mag); % align output, see note below
figure;
stem(F,yfft_mag);
xlabel('Frequency (Fs^{-1})');
%fftshift() cuts the second (complex) half of the spectrum and pastes it backwards at the beginning, so that our frequency
%axis is now correct; it is in units of 1/Fs, so 0.1 corresponds to the 2Hz we put in

keyboard;


%% FFT output - magnitude and phase spectra

% Complete periods of signal not present:
tvec = t0:1/Fs:t1;
y = sin(2*pi*f*tvec);
yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
figure;
stem(yfft_mag);
% No frequency bin that is true freq of signal.

% Control nuber of points in fft
% Zero padding allows a longe harmonic series with smaller fundamental frequency
tvec = t0:1/Fs:t1;
nPoints = [length(tvec) 64 256 1024];
figure; hold on;
for iP = 1:length(nPoints) % repeat fft with different numbers of points
    nP = nPoints(iP);
    subplot(2,2,iP);
    y = sin(2*pi*f*tvec);
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    title(sprintf('%d point FFT',nP));
    xlabel('Frequency (Fs^{-1})');
end
% Zero padding original signal will also have similar effect

keyboard;

%% Spectral leakage
% ------------------

tvec = t0:1/Fs:t1-(1/Fs);
nRepeats = [1 2 4 8];
nP = 1024;
figure; hold on;

for iP = 1:length(nRepeats)
    subplot(2,2,iP);
    y = sin(2*pi*f*tvec);
    y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
    title(sprintf('%d repeats',nRepeats(iP)));
    xlabel('Frequency (Fs^{-1})');
end


keyboard;

%% WINDOWING
% ----------
% Thereore, use windowing to prevent abrupt changes in signals

% CONVOLUTION THEOREM! - HALT!!



nP = 25;
nPFFT = 1024;
windows = {'rectwin','triang','hamming','hanning','blackman'};
cols = 'rgbcmyk';
figure; hold on;
for iW = 1:length(windows)
    eval(cat(2,'wn = ',windows{iW},'(nP);')); % make sure you understand this
    wn = wn./sum(wn);
    subplot(211); % plot the window
    plot(wn,cols(iW),'LineWidth',2); hold on;
    subplot(212);
    yfft = fft(wn,nPFFT);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
    F = [-nPFFT/2:nPFFT/2-1]./nPFFT;
    yfft_mag = fftshift(yfft_mag);
    h(iW) = plot(F,yfft_mag,cols(iW),'LineWidth',2); hold on;
end
xlabel('Frequency (Fs^{-1})');
legend(h,windows);


%% Spectral Estimation

% periodogram // get Fourier coefficients by aplying DFT to entire data set
% Noisy and biased

[Pxx,F] = periodogram(y,[],nP,Fs);
figure;
plot(F,Pxx); xlabel('Frequency (Hz)');
hold on;
[Pxx,F] = periodogram(y,hanning(length(y)),nP,Fs);
plot(F,Pxx,'r');

keyboard;

% More robust estimates by estimating spectrums for segments - but tradeoff
% in frequency resolution
% Bartletts method/ Welch's method 

t0 = 0; t1 = 1;
f = 2;
nRepeats = 4;
tvec = t0:1/Fs:t1-(1/Fs);
nP = 1024;
y = sin(2*pi*f*tvec);
y = repmat(y,[1 nRepeats]);
[Pxx,F] = periodogram(y,rectwin(length(y)),nP,Fs);
figure ; hold on
plot(F,Pxx);
hold on;
wSize = 40;
[Pxx,F] = pwelch(y,rectwin(wSize),wSize/2,nP,Fs);
plot(F,Pxx,'r'); xlabel('Frequency (Hz)');


keyboard; 

%% Unevenly sampled data
% ----------------------

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; f = 2;
nP = 1024;
gaps = [5 10 15]; % idx of samples to be removed
tvec = t0:1/Fs:t1;%-(1/Fs);
y = sin(2*pi*f*tvec);
figure; hold on;
subplot(211)
plot(tvec,y,'k*'); hold on;
yfft = fft(y,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
subplot(212);
plot(F,yfft_mag,'k-x'); hold on;
xlabel('Frequency (Fs^{-1})');

% signal with gaps

y2 = y;
y2(gaps) = []; tvec(gaps) = []; % remove
figure; hold on;
subplot(211);
plot(tvec,y2,'bo'); hold on;
yfft = fft(y2,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
subplot(212);
plot(F,yfft_mag,'b-x');
legend('y','y2 (with gaps)')


