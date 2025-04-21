
% Sampling - Aliasing
% --------------------
% Also take a look at the Sampling Handouts!


set(0,'DefaultAxesFontSize',18)

% Sampling a 10Hz sinusoidal signal at 1000 Hz

Fs1 = 1000; % Fs is the conventional variable name for sampling freq
F1 = 10; twin = [0 1]; % use a 1-second time window (from 0 to 1s)

tvec1 = twin(1):1/Fs1:twin(2); % timebase for signal
signal1 = sin(2*pi*F1*tvec1);

% Sample at 12 Hz
Fs2 = 12;
tvec2 = twin(1):1/Fs2:twin(2);
signal2 = interp1(tvec1,signal1,tvec2,'nearest');
% interp1 - inout signal specified by tvec1 and signal 1 - return for tvec2, "nearest" values in tvec1

figure; hold on;
plot(tvec1,signal1);
hold on;
plot(tvec2,signal2,'.g','MarkerSize',20);
plot(tvec1,-sin(2*pi*2*tvec1),'r--','LineWidth',2);
xlabel('time (s)'); ylabel('y');


keyboard;

% Aliasing (sub-sampling/ decimating )
% --------------------------------------
% Nyquist - sampling freq should be at-least twice of signal frequency
% Use anti-alasing filter

Fs1 = 1200;
F1 = 3; F2 = 10;  
twin = [0 1];
tvec1 = twin(1):1/Fs1:twin(2);

signal1a = sin(2*pi*F1*tvec1); signal1b = 0.5*sin(2*pi*F2*tvec1);
signal1 = signal1a + signal1b;

% Sub-sample at 12 Hz
figure; hold on;
dt = 100;
tvec2 = tvec1(1:dt:end);
signal2 = signal1(1:dt:end); % sample at 12 Hz - every 100th sample
subplot(131)
plot(tvec1,signal1);
hold on;
stem(tvec2,signal2,'.g','MarkerSize',20);
title('without anti-aliasing filter'); % spurious 10Hz signal is aliased

% sample at 12 Hz with different method
tvec1d = decimate(tvec1, dt);
signal2d = decimate(signal1,dt);      % decimate filters the data
subplot(132)
plot(tvec1,signal1a,'b--');
hold on;
stem(tvec1d,signal2d,'.g','MarkerSize',20);
xlabel('time (s)'); ylabel('y');
title('with anti-aliasing filter');
subplot(133)
plot(tvec2,signal2-signal2d,'r','LineWidth',2);
title('difference');


keyboard;

% reconstruct signal from sampled data
% ------------------------------------

figure; hold on;
fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
ax1 = subplot(211);
stem(tvec,y); title('original');

subsample_factor = 4;
tvec2 = tvec(1:subsample_factor:end); % take every 4th sample
y2 = y(1:subsample_factor:end);
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');

xl = [1 1.04];
linkaxes([ax1, ax2], 'x');
set(ax1,'XLim',xl); %

hold on;
y_interp = interp1(tvec2,y2,tvec,'linear');
p1 = plot(tvec,y_interp,'b');
y_interp2 = interp1(tvec2,y2,tvec,'spline');
p2 = plot(tvec,y_interp2,'g');
legend([p1 p2],{'linear','spline'},'Location','Northeast'); legend boxoff


keyboard;











