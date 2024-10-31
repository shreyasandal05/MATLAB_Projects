clear
close all
clc

Ts=1e-3;% 1 ms
sps = 16;  % Number of samples per symbol (oversampling factor)
T_sample=Ts/sps; % sps samples in each symbol duration
F_sample=1/T_sample;
fc = 100;
M = 16;   % 16-QAM
k = log2(M); % Bits per symbol
SNR = 20; % Signal-to-Noise Ratio in dB
numBits = 1000;
num_symbols = numBits / log2(M);

% Generate random bits
dataBits = randi([0 1], 1, numBits);


%% Grouping and Mapping

% Convert bits to symbols
dataSymbols = bi2de(reshape(dataBits, k, []).', 'left-msb');

% Map symbols to IQ points
% 16-QAM Constellation (example mapping)
constellation = [
    -3 -3;  -3 -1;  -3 +1;  -3 +3;
    -1 -3;  -1 -1;  -1 +1;  -1 +3;
    +1 -3;  +1 -1;  +1 +1;  +1 +3;
    +3 -3;  +3 -1;  +3 +1;  +3 +3
];

% Generate IQ samples based on the symbols
IQ_Samples = zeros(1, num_symbols * sps);
for i = 1:num_symbols
    symbol = dataSymbols(i) + 1;  % Convert to 1-based index
    IQ_Samples((i-1)*sps + 1:i*sps) = ...
        constellation(symbol, 1) + 1j * constellation(symbol, 2);
end

% Time vector for plotting
t = 0:Ts:(length(IQ_Samples)-1)*Ts;

% Plot the IQ samples
figure(1);
subplot(2, 1, 1);
plot(t, real(IQ_Samples)); hold on
title('Real Part of IQ Samples');
grid on

subplot(2, 1, 2);
plot(t, imag(IQ_Samples)); hold on
title('Imaginary Part of IQ Samples');
grid on

% Plot the constellation
figure(2);
scatter(constellation(:, 1), constellation(:, 2), 'filled');
title('16-QAM Constellation Diagram');
axis equal;
grid on

%% Pulse shapping

% Pulse Shaping using a Root Raised Cosine Filter
rolloff = 0.25;               % Rolloff factor
span = 6;                     % Span of the filter in symbols
r = rcosdesign(rolloff, span, sps);  % Generate the filter

% Convolve IQ samples with the pulse shaping filter
shapedSignal = conv(IQ_Samples, r, 'same');

%% %% PASSBAND CONVERSION

% Passband Modulation
t = (0:length(shapedSignal)-1) / F_sample;  % Time vector
I_passbandSignal = real(shapedSignal) .* cos(2 * pi * fc * t);
Q_passbandSignal = -imag(shapedSignal) .* sin(2 * pi * fc * t);

%% passband IQ samples
figure(3)
subplot(2,1,1)
plot(t,real(shapedSignal), 'r-'); hold on
plot(t,-real(shapedSignal), 'r-'); hold on
plot(t,I_passbandSignal, 'b.-'); hold on
title("I Samples in passband")
subplot(2,1,2)
plot(t,imag(shapedSignal), 'r-'); hold on
plot(t,-imag(shapedSignal), 'r-'); hold on
plot(t,Q_passbandSignal, 'b.-'); hold on
title("Q Samples in passband")


%% Noise Adding

% Add AWGN noise

noisySignal = awgn(shapedSignal, SNR, 'measured');

% Demodulation (assuming coherent detection)
% Sample the received signal
receivedSamples = noisySignal(1:sps:end);  % Take one sample per symbol



% Plotting
figure(4);
subplot(1, 1, 1);
plot(real(shapedSignal));
title('Pulse Shaped Signal');
xlabel('Samples');
ylabel('Amplitude');

figure(5)
subplot(1, 1, 1);
plot(real(noisySignal));
title('Received Noisy Signal');
xlabel('Samples');
ylabel('Amplitude');

figure(6)
subplot(1, 1, 1);
scatter(constellation(:, 1), constellation(:, 2), 'filled');
hold on;
scatter(real(receivedSamples), imag(receivedSamples), 'r');  % Received symbols
title('Constellation Diagram');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal;
grid on;

%% % --- Passband to Baseband Conversion ---

% Demodulate the Passband Signal
x_I_baseband_down = I_passbandSignal .* cos(2 * pi * fc * t);  % I component
x_Q_baseband_down = Q_passbandSignal .* sin(2 * pi * fc * t);  % Q component


% Design a low-pass filter to recover the I and Q components
Rb = 1/Ts;
rolloff = 0.25;
wc = Rb/2*(1+rolloff);
b = fir1(50,2*pi*wc/F_sample);
y_I_shaped = filter(b,1,x_I_baseband_down);
y_Q_shaped = filter(b,1,x_Q_baseband_down);

% Apply the low-pass filter to both components
 x_I_baseband_down = filter(b, 1, y_I_shaped );
x_Q_baseband_down = filter(b, 1, y_Q_shaped );

% Downsample to recover the baseband signal (taking one sample per symbol)
baseband_I = x_I_baseband_down(1:sps:end);
baseband_Q = x_Q_baseband_down(1:sps:end);

% Combine I and Q components into complex baseband signal
basebandSignal = baseband_I + 1j * baseband_Q;

% Map back to symbols
demodulatedSymbols = zeros(1, num_symbols);
for i = 1:num_symbols
    % Find closest constellation points
    [~, demodulatedSymbols(i)] = min(abs(basebandSignal(i) - constellation(:, 1) + 1j * constellation(:, 2)));
end

% Convert symbols back to bits
demodulatedBits = de2bi(demodulatedSymbols - 1, log2(M), 'left-msb');  % Convert back to bits
demodulatedBits = demodulatedBits(:)';  % Reshape to a row vector

% Plotting the Signals
figure(7);
subplot(2, 1, 1);
plot(t, x_I_baseband_down , 'b'); % Plot I component
hold on;
plot(t, x_Q_baseband_down , 'r--'); % Plot Q component
title('Demodulated I and Q Components');
xlabel('Time (s)');
ylabel('Amplitude');
legend('I Component', 'Q Component');
grid on;

subplot(2, 1, 2);
scatter(constellation(:, 1), constellation(:, 2), 'filled');
hold on;
scatter(real(basebandSignal), imag(basebandSignal), 'r');  % Demodulated symbols
title('Constellation Diagram after Demodulation');
xlabel('In-Phase');
ylabel('Quadrature');
axis equal;
grid on;

% Calculate Bit Error Rate (BER)
ber = sum(dataBits ~= demodulatedBits) / numBits;  % Bit Error Rate

% Display Results
disp(['Bit Error Rate (BER): ', num2str(ber)]);
%% SNR

filtlen = 10;      % Filter length in symbols


rrcFilter = rcosdesign(rolloff,filtlen,sps);
impz(rrcFilter)

rng default;                     % Default random number generator
dataIn = randi([0 1],numBits,1); % Generate vector of binary data
dataSymbolsIn = bit2int(dataIn,k);
dataMod = qammod(dataSymbolsIn,M);

txFiltSignal = upfirdn(dataMod,rrcFilter,sps,1);

EbNo = 10;
snr = convertSNR(EbNo,'ebno', ...
    samplespersymbol=sps, ...
    bitspersymbol=k);
rxSignal = awgn(txFiltSignal,snr,'measured');

rxFiltSignal = ...
    upfirdn(rxSignal,rrcFilter,1,sps);       % Downsample and filter
rxFiltSignal = ...
    rxFiltSignal(filtlen + 1:end - filtlen); % Account for delay

dataSymbolsOut = qamdemod(rxFiltSignal,M);
dataOut = int2bit(dataSymbolsOut,k);

[numErrors,ber] = biterr(dataIn,dataOut);
fprintf(['\nFor an EbNo setting of %3.1f dB, ' ...
    'the bit error rate is %5.2e, based on %d errors.\n'], ...
    EbNo,ber,numErrors)

EbNo = 10;
snr = convertSNR(EbNo,'ebno', ...
    samplespersymbol=sps, ...
    bitspersymbol=k);
rxSignal = awgn(txFiltSignal,snr,'measured');
rxFiltSignal = ...
    upfirdn(rxSignal,rrcFilter,1,sps);       % Downsample and filter
rxFiltSignal = ...
    rxFiltSignal(filtlen + 1:end - filtlen); % Account for delay

eyediagram(txFiltSignal(1:200),sps*2);
eyediagram(rxSignal(1:200),sps*2);
eyediagram(rxFiltSignal(1:200),2);


scatplot = scatterplot(sqrt(sps)*...
    rxSignal(1:sps*5e3),...
    sps,0);
hold on;
scatterplot(rxFiltSignal(1:5e3),1,0,'bx',scatplot);
title('Received Signal, Before and After Filtering');
legend('Before Filtering','After Filtering');
axis([-5 5 -5 5]); % Set axis ranges
hold off;