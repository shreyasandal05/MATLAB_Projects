% Parameters
Fs = 1000;                 % Sampling frequency in Hz
T = 1;                     % Duration in seconds
t = 0:1/Fs:T-1/Fs;         % Time vector
data = randi([0 1], 1, 20); % Random binary data sequence
n = length(data);          % Number of bits
bitDuration = T/n;         % Duration of one bit

% Unipolar NRZ Signal
nrzSignal = [];
for i = 1:n
    if data(i) == 1
        nrzSignal = [nrzSignal ones(1, Fs*bitDuration)];
    else
        nrzSignal = [nrzSignal zeros(1, Fs*bitDuration)];
    end
end
t_nrz = 0:1/Fs:(length(nrzSignal)-1)/Fs;

% Unipolar RZ Signal
rzSignal = [];
for i = 1:n
    if data(i) == 1
        rzSignal = [rzSignal ones(1, Fs*bitDuration/2) zeros(1, Fs*bitDuration/2)];
    else
        rzSignal = [rzSignal zeros(1, Fs*bitDuration)];
    end
end
t_rz = 0:1/Fs:(length(rzSignal)-1)/Fs;

% Compute and Plot the Spectrum
nfft = 2048; % FFT size
f = Fs*(0:(nfft/2))/nfft; % Frequency vector

% Spectrum of Unipolar NRZ
nrzSpectrum = abs(fft(nrzSignal, nfft));
nrzSpectrum = nrzSpectrum(1:nfft/2+1);

% Spectrum of Unipolar RZ
rzSpectrum = abs(fft(rzSignal, nfft));
rzSpectrum = rzSpectrum(1:nfft/2+1);

% Power Spectral Density using Welch's method
[psdNrz, f_psdNrz] = pwelch(nrzSignal, [], [], [], Fs);
[psdRz, f_psdRz] = pwelch(rzSignal, [], [], [], Fs);

% Plot Signals and Spectra
figure;

% Plot Unipolar NRZ Signal
subplot(2,2,1);
plot(t_nrz, nrzSignal, 'LineWidth', 1.5);
title('Unipolar NRZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -0.1 1.1]);

% Plot Unipolar RZ Signal
subplot(2,2,2);
plot(t_rz, rzSignal, 'LineWidth', 1.5);
title('Unipolar RZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -0.1 1.1]);

% Plot Spectrum of Unipolar NRZ
subplot(2,2,3);
plot(f, 20*log10(nrzSpectrum));
title('Spectrum of Unipolar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Plot Spectrum of Unipolar RZ
subplot(2,2,4);
plot(f, 20*log10(rzSpectrum));
title('Spectrum of Unipolar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Create a new figure for Power Spectral Density
figure;

% Plot Power Spectral Density of Unipolar NRZ
subplot(2,1,1);
plot(f_psdNrz, 10*log10(psdNrz));
title('Power Spectral Density of Unipolar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Plot Power Spectral Density of Unipolar RZ
subplot(2,1,2);
plot(f_psdRz, 10*log10(psdRz));
title('Power Spectral Density of Unipolar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Adjust plot settings
sgtitle('Unipolar NRZ and RZ Signals with Spectra and Power Spectral Density');
