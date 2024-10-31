% Parameters
numBits = 64;                  % Total number of bits (must be a multiple of 4 for 16-QAM)
M = 16;                        % Modulation order (16-QAM)
numSymbols = numBits / log2(M);  % Number of symbols
sps = 4;                       % Samples per symbol
rolloff = 0.25;               % Rolloff factor for pulse shaping
span = 6;                     % Span of the filter in symbols
fs = 1000;                    % Sampling frequency for baseband signal
EbN0_dB = 10;                 % Energy per bit to noise power spectral density ratio in dB

% Generate random bits
dataBits = randi([0 1], 1, numBits);

% Convert bits to symbols
dataSymbols = bi2de(reshape(dataBits, log2(M), []).', 'left-msb');

% 16-QAM Constellation Mapping
constellation = [
    -3 -3;  -3 -1;  -3 +1;  -3 +3;
    -1 -3;  -1 -1;  -1 +1;  -1 +3;
    +1 -3;  +1 -1;  +1 +1;  +1 +3;
    +3 -3;  +3 -1;  +3 +1;  +3 +3
];

% Generate IQ Samples
IQ_Samples = zeros(1, numSymbols * sps);
for i = 1:numSymbols
    symbol = dataSymbols(i) + 1;  % Convert to 1-based index
    IQ_Samples((i-1)*sps + 1:i*sps) = ...
        constellation(symbol, 1) + 1j * constellation(symbol, 2);
end

% Pulse Shaping using a Root Raised Cosine Filter
r = rcosdesign(rolloff, span, sps);  % Generate the RRC filter
shapedSignal = conv(IQ_Samples, r, 'same');

% Add AWGN Noise
snr = EbN0_dB + 10*log10(log2(M)) - 10*log10(sps); % Calculate SNR in dB
noisySignal = awgn(shapedSignal, snr, 'measured');

% --- Demodulation Process ---
% Downsample to recover the baseband signal
downsampledSignal = noisySignal(1:sps:end);

% Map back to symbols
demodulatedSymbols = zeros(1, numSymbols);
for i = 1:numSymbols
    % Find closest constellation points
    [~, demodulatedSymbols(i)] = min(abs(downsampledSignal(i) - constellation(:, 1) + 1j * constellation(:, 2)));
end

% Convert symbols back to bits
demodulatedBits = de2bi(demodulatedSymbols - 1, log2(M), 'left-msb');  % Convert back to bits
demodulatedBits = demodulatedBits(:)';  % Reshape to a row vector

% Calculate Bit Error Rate (BER)
ber = sum(dataBits ~= demodulatedBits) / numBits;  % Bit Error Rate

% Display Results
disp(['Bit Error Rate (BER): ', num2str(ber)]);

% Plotting the Signals
figure(1);
plot(real(shapedSignal));
title('Transmitted Baseband Signal');
xlabel('Samples');
ylabel('Amplitude');
grid on;

figure(2);
plot(real(noisySignal));
title('Received Noisy Baseband Signal');
xlabel('Samples');
ylabel('Amplitude');
grid on;

figure(3);
scatter(constellation(:, 1), constellation(:, 2), 'filled');
hold on;
scatter(real(downsampledSignal), imag(downsampledSignal), 'r');  % Demodulated symbols
title('Constellation Diagram after Demodulation');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal;
grid on;