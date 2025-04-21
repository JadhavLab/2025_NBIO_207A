%% Spectrograms
% -----------------------------------------


% data!
%% construct and plot the spectrogram
[S,F,T,P] = spectrogram(cscR.data,hanning(512),256,1:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');

% we specified 512-sample windows with a 256-sample overlap, we expect roughly 8 estimates per
%second. This is because our Fs is 2000, and therefore 8 ?steps? of our moving window, moving at 256 samples 
% per step (the window size minus the overlap equals the step size) fit into one second

% the time of the first bin returned by spectrogram() is the center of the first window for which there is enough data.

% Spectrogram parameters
% --------------------

%Change the code above so that spectrogram() uses an overlap of 384 samples (instead of 256) and returns power with 0.25 Hz resolution
%(instead of 1).

% Change the window used to a rectangular one (rectwin). What happens to the spectrogram?

% compare side by side (using subplot()) the spectrogram obtained with 512- and 1024-sample windows. For a
% fair comparison, set the overlap such that you end up with a comparable number of time bins.
