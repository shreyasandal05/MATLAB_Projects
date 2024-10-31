clear
close all
clc

Ts=1e-3;% 1 ms
sps = 16;
T_sample=Ts/sps;          % sps samples in each symbol duration
F_sample=1/T_sample; 

num_symbols = 100;
bits = randi([0, 1],1, num_symbols*2); %Our data to be transmitted, 1's and 0's
bits_I = bits(1:2:end);
bits_Q = bits(2:2:end);

%% BASEBAND PROCESSING
%% In-Phase
x_I = [];
t=(0:sps*num_symbols-1)*T_sample;
for bit=bits_I
    pulse = zeros(1,sps);
    pulse(1) = bit*2-1; %set the first value to either a 1 or -1
    x_I = [x_I, pulse]; %add the 8 samples to the signal
end

figure(1)
stem(t, x_I, 'b-o'); hold on
grid on 

%% Quadrature
x_Q = [];
t=(0:sps*num_symbols-1)*T_sample;
for bit=bits_Q
    pulse = zeros(1,sps);
    pulse(1) = bit*2-1; %set the first value to either a 1 or -1
    x_Q = [x_Q, pulse]; %add the 8 samples to the signal
end

figure(1)
stem(t, x_Q, 'r-*')
title("BB - IQ Samples")
grid on

% Create our root-raised-cosine filter
filtlen = 10;      % Filter length in symbols
beta = 0.25;    % Filter beta factor
h = rcosdesign(beta,filtlen,sps);
num_taps=length(h);
figure(2)
plot(h, '.')
title("RRC Samples")
grid on

%% Pulse Shaping 

% I-channel
x_shaped_I = conv(x_I, h);
% Q-channel
x_shaped_Q = conv(x_Q, h);
t=(0:length(x_shaped_I)-1)*T_sample;

figure(3)
plot(t,x_shaped_I, 'b.-',t,x_shaped_Q, 'r.-'); hold on;grid on
title("Pulse shaped IQ Samples")
legend('I-Chan','Q-Chan')

%% PASSBAND CONVERSION

Fc=5000; % Carrier frequency
x_I_passband=x_shaped_I.*cos(2*pi*Fc*t)*sqrt(2);
x_Q_passband=-x_shaped_Q.*sin(2*pi*Fc*t)*sqrt(2);

%% passband IQ samples
figure(4)
subplot(2,1,1)
plot(t,x_shaped_I, 'r-'); hold on
plot(t,-x_shaped_I, 'r-'); hold on
plot(t,x_I_passband, 'b.-'); hold on
title("IQ Samples in passband")
subplot(2,1,2)
plot(t,x_shaped_Q, 'r-'); hold on
plot(t,-x_shaped_Q, 'r-'); hold on
plot(t,x_Q_passband, 'b.-'); hold on
title("IQ Samples in passband")


%% Power spectrum in baseband
figure(5)
subplot(2,2,1)
[pxx,f] = pspectrum(x_shaped_I, F_sample,'Leakage',1,'TwoSided',true);
plot(f,pxx);grid on
title("Power spectrum in  I-baseband")
subplot(2,2,3)
[pxx,f] = pspectrum(x_shaped_Q, F_sample,'Leakage',1,'TwoSided',true);
plot(f,pxx);grid on
title("Power spectrum in  Q-baseband")

%% Power spectrum in passband
figure(5)
subplot(2,2,2)
[pxx,f] = pspectrum(x_I_passband, F_sample,'Leakage',1,'TwoSided',true);
plot(f,pxx);grid on
title("Power spectrum I in passband")
subplot(2,2,4)
[pxx,f] = pspectrum(x_Q_passband, F_sample,'Leakage',1,'TwoSided',true);
plot(f,pxx);grid on
title("Power spectrum Q in passband")

%% Transmit Signal
x_t_passband=x_I_passband+x_Q_passband;
figure(6)
subplot(2,1,1)
plot(t,x_t_passband, 'b-'); hold on
title('trasmit signal s(t)')
figure (6)
subplot(2,1,2)
[pxx,f] = pspectrum(x_t_passband, F_sample,'Leakage',1,'TwoSided',true);
plot(f,pxx);grid on
title('Spectrum of trasmit signal S(f)')

%% Add AWGN CHANNEL
noise=0.03*randn(size(x_t_passband));
x_noisy=x_t_passband+noise;
figure(6)
subplot(2,1,1)
plot(t,x_noisy, 'r-'); hold on
title('received signal s(t)+n(t)')

%% Receiver processing
% PASSBAND TO BASEBAND

x_I_baseband_down = x_noisy.*cos(2*pi*Fc*t)*sqrt(2);
x_Q_baseband_down = -x_noisy.*sin(2*pi*Fc*t)*sqrt(2);
Rb = 1/Ts;
wc = Rb/2*(1+beta);

%% Low-pass filter design
b = fir1(50,2*pi*wc/F_sample);
y_I_shaped = filter(b,1,x_I_baseband_down);
y_Q_shaped = filter(b,1,x_Q_baseband_down);

figure(7)
subplot(2,1,1)
plot(t,x_shaped_I, 'r-'); hold on;grid on;
plot(t,y_I_shaped, 'b-'); hold on
title('down converted I bb')
figure(7)
subplot(2,1,2)
plot(t,x_shaped_Q, 'r-'); hold on;grid on;
plot(t,y_Q_shaped, 'b-'); hold on
title('down converted Q bb')
%% MATCHED FILTER
y_I_received = conv(y_I_shaped, h);
y_Q_received = conv(y_Q_shaped, h);

%plot matched filter output
figure(8)
subplot(2,1,1)
t=(0:length(y_I_received)-1)*T_sample;
plot(t,y_I_received, 'r.-'); hold on;grid on
title('matched filter I chann output')

subplot(2,1,2)
plot(t,y_Q_received, 'r.-'); hold on;grid on
title('matched filter Q chann output')

% Sample the signal at Ts, consider the delay in filter
figure(8)
subplot(2,1,1)
nn=(0:num_symbols-1)*sps+(num_taps-1)+(length(b)-1)/2+1;
y_I=y_I_received(nn)
stem((nn-1)*T_sample,y_I,'go',"filled");

subplot(2,1,2)
y_Q=y_Q_received(nn)
stem((nn-1)*T_sample,y_Q,'go',"filled");
scatterplot(y_I+1j*y_Q)