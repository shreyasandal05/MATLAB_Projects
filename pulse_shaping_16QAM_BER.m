clear

M = 16;            % Modulation order
k = log2(M);       % Bits per symbol
numBits = k*7.5e4; % Bits to process
sps = 4;           % Samples per symbol (oversampling factor)

filtlen = 10;      % Filter length in symbols
rolloff = 0.25;    % Filter rolloff factor

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

EbNo = 20;
snr = convertSNR(EbNo,'ebno', ...
    samplespersymbol=sps, ...
    bitspersymbol=k);
rxSignal = awgn(txFiltSignal,snr,'measured');
rxFiltSignal = ...
    upfirdn(rxSignal,rrcFilter,1,sps);       % Downsample and filter
rxFiltSignal = ...
    rxFiltSignal(filtlen + 1:end - filtlen); % Account for delay

eyediagram(txFiltSignal(1:2000),sps*2);
eyediagram(rxSignal(1:2000),sps*2);
eyediagram(rxFiltSignal(1:2000),2);


scatplot = scatterplot(sqrt(sps)*...
    rxSignal(1:sps*5e3),...
    sps,0);
hold on;
scatterplot(rxFiltSignal(1:5e3),1,0,'bx',scatplot);
title('Received Signal, Before and After Filtering');
legend('Before Filtering','After Filtering');
axis([-5 5 -5 5]); % Set axis ranges
hold off;