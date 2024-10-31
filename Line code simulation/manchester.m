% Parameters
Fs = 1000;                 % Sampling frequency in Hz
T = 1;                     % Duration in seconds
data = randi([0 1], 1, 20); % Random binary data sequence
n = length(data);          % Number of bits
bitDuration = T/n;         % Duration of one bit
t = 0:1/Fs:T-1/Fs;         % Time vector for full duration

% Manchester Encoding
manchesterSignal = [];
for i = 1:n
    if data(i) == 1
        % For binary 1: high-to-low transition in the middle
        manchesterSignal = [manchesterSignal ones(1, Fs*bitDuration/2) -ones(1, Fs*bitDuration/2)];
    else
        % For binary 0: low-to-high transition in the middle
        manchesterSignal = [manchesterSignal -ones(1, Fs*bitDuration/2) ones(1, Fs*bitDuration/2)];
    end
end
t_manchester = 0:1/Fs:(length(manchesterSignal)-1)/Fs;

% Compute Spectrum
nfft = 2048; % FFT size
f = Fs*(0:(nfft/2))/nfft; % Frequency vector

% Spectrum of Manchester Signal
manchesterSpectrum = abs(fft(manchesterSignal, nfft));
manchesterSpectrum = manchesterSpectrum(1:nfft/2+1);

% Power Spectral Density using Welch's method
[psdManchester, f_psdManchester] = pwelch(manchesterSignal, [], [], [], Fs);

% Plot Signals and Spectra
figure;

% Plot Manchester Encoded Signal
subplot(2,2,1);
plot(t_manchester, manchesterSignal, 'LineWidth', 1.5);
title('Manchester Encoded Signal');
xlabel('Time (s)');
ylabel('Amplitude');
axis([0 T -1.5 1.5]);

% Plot Spectrum of Manchester Signal
subplot(2,2,2);
plot(f, 20*log10(manchesterSpectrum));
title('Spectrum of Manchester Encoded Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

% Plot Power Spectral Density of Manchester Signal
subplot(2,1,2);
plot(f_psdManchester, 10*log10(psdManchester));
title('Power Spectral Density of Manchester Encoded Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Adjust plot settings
sgtitle('Manchester Encoding: Signal, Spectrum, and Power Spectral Density');
