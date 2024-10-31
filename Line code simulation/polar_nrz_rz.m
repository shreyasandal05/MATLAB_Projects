% Parameters
Fs = 1000;                 % Sampling frequency in Hz
T = 1;                     % Duration in seconds
t = 0:1/Fs:T-1/Fs;         % Time vector
data = randi([0 1], 1, 20); % Random binary data sequence
n = length(data);          % Number of bits
bitDuration = T/n;         % Duration of one bit

% Polar NRZ Signal
polarNrzSignal = [];
for i = 1:n
    if data(i) == 1
        polarNrzSignal = [polarNrzSignal ones(1, Fs*bitDuration)];
    else
        polarNrzSignal = [polarNrzSignal -ones(1, Fs*bitDuration)];
    end
end
t_polar_nrz = 0:1/Fs:(length(polarNrzSignal)-1)/Fs;

% Polar RZ Signal
polarRzSignal = [];
for i = 1:n
    if data(i) == 1
        polarRzSignal = [polarRzSignal ones(1, Fs*bitDuration/2) -ones(1, Fs*bitDuration/2)];
    else
        polarRzSignal = [polarRzSignal -zeros(1, Fs*bitDuration)];
    end
end
t_polar_rz = 0:1/Fs:(length(polarRzSignal)-1)/Fs;

% Compute and Plot the Spectrum
nfft = 2048; % FFT size
f = Fs*(0:(nfft/2))/nfft; % Frequency vector

% Spectrum of Polar NRZ
polarNrzSpectrum = abs(fft(polarNrzSignal, nfft));
polarNrzSpectrum = polarNrzSpectrum(1:nfft/2+1);

% Spectrum of Polar RZ
polarRzSpectrum = abs(fft(polarRzSignal, nfft));
polarRzSpectrum = polarRzSpectrum(1:nfft/2+1);

% Power Spectral Density using Welch's method
[psdNrz, f_psdNrz] = pwelch(polarNrzSignal, [], [], [], Fs);
[psdRz, f_psdRz] = pwelch(polarRzSignal, [], [], [], Fs);

% Plot Signals and Spectra
figure;

% Plot Polar NRZ Signal
subplot(2,2,1);
plot(t_polar_nrz, polarNrzSignal, 'LineWidth', 1.5);
title('Polar NRZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -1.5 1.5]);

% Plot Polar RZ Signal
subplot(2,2,2);
plot(t_polar_rz, polarRzSignal, 'LineWidth', 1.5);
title('Polar RZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -1.5 1.5]);

% Plot Spectrum of Polar NRZ
subplot(2,2,3);
plot(f, 20*log10(polarNrzSpectrum));
title('Spectrum of Polar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Plot Spectrum of Polar RZ
subplot(2,2,4);
plot(f, 20*log10(polarRzSpectrum));
title('Spectrum of Polar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Create a new figure for Power Spectral Density
figure;

% Plot Power Spectral Density of Polar NRZ
subplot(2,1,1);
plot(f_psdNrz, 10*log10(psdNrz));
title('Power Spectral Density of Polar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Plot Power Spectral Density of Polar RZ
subplot(2,1,2);
plot(f_psdRz, 10*log10(psdRz));
title('Power Spectral Density of Polar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Adjust plot settings
sgtitle('Polar NRZ and RZ Signals with Spectra and Power Spectral Density');
