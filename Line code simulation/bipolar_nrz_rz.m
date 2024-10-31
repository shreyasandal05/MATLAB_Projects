% Parameters
Fs = 1000;                 % Sampling frequency in Hz
T = 1;                     % Duration in seconds
t = 0:1/Fs:T-1/Fs;         % Time vector
data = randi([0 1], 1, 20); % Random binary data sequence
n = length(data);          % Number of bits
bitDuration = T/n;         % Duration of one bit

% Bipolar NRZ Signal
bipolarNrzSignal = [];
for i = 1:n
    if data(i) == 1
        bipolarNrzSignal = [bipolarNrzSignal ones(1, Fs*bitDuration) -ones(1, Fs*bitDuration)];
    else
        bipolarNrzSignal = [bipolarNrzSignal zeros(1, Fs*bitDuration) zeros(1, Fs*bitDuration)];
    end
end
t_bipolar_nrz = 0:1/Fs:(length(bipolarNrzSignal)-1)/Fs;

% Bipolar RZ Signal
bipolarRzSignal = [];
for i = 1:n
    if data(i) == 1
        bipolarRzSignal = [bipolarRzSignal ones(1, Fs*bitDuration/2) -ones(1, Fs*bitDuration/2)];
    else
        bipolarRzSignal = [bipolarRzSignal zeros(1, Fs*bitDuration)];
    end
end
t_bipolar_rz = 0:1/Fs:(length(bipolarRzSignal)-1)/Fs;

% Compute and Plot the Spectrum
nfft = 2048; % FFT size
f = Fs*(0:(nfft/2))/nfft; % Frequency vector

% Spectrum of Bipolar NRZ
bipolarNrzSpectrum = abs(fft(bipolarNrzSignal, nfft));
bipolarNrzSpectrum = bipolarNrzSpectrum(1:nfft/2+1);

% Spectrum of Bipolar RZ
bipolarRzSpectrum = abs(fft(bipolarRzSignal, nfft));
bipolarRzSpectrum = bipolarRzSpectrum(1:nfft/2+1);

% Power Spectral Density using Welch's method
[psdNrz, f_psdNrz] = pwelch(bipolarNrzSignal, [], [], [], Fs);
[psdRz, f_psdRz] = pwelch(bipolarRzSignal, [], [], [], Fs);

% Plot Signals and Spectra
figure;

% Plot Bipolar NRZ Signal
subplot(2,2,1);
plot(t_bipolar_nrz, bipolarNrzSignal, 'LineWidth', 1.5);
title('Bipolar NRZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -1.5 1.5]);

% Plot Bipolar RZ Signal
subplot(2,2,2);
plot(t_bipolar_rz, bipolarRzSignal, 'LineWidth', 1.5);
title('Bipolar RZ Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -1.5 1.5]);

% Plot Spectrum of Bipolar NRZ
subplot(2,2,3);
plot(f, 20*log10(bipolarNrzSpectrum));
title('Spectrum of Bipolar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Plot Spectrum of Bipolar RZ
subplot(2,2,4);
plot(f, 20*log10(bipolarRzSpectrum));
title('Spectrum of Bipolar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Create a new figure for Power Spectral Density
figure;

% Plot Power Spectral Density of Bipolar NRZ
subplot(2,1,1);
plot(f_psdNrz, 10*log10(psdNrz));
title('Power Spectral Density of Bipolar NRZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Plot Power Spectral Density of Bipolar RZ
subplot(2,1,2);
plot(f_psdRz, 10*log10(psdRz));
title('Power Spectral Density of Bipolar RZ Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Adjust plot settings
sgtitle('Bipolar NRZ and RZ Signals with Spectra and Power Spectral Density');
