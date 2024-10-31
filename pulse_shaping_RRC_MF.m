clear 

Ts=1e-3;% 1 ms
sps = 8;
T_sample=Ts/sps;          % sps samples in each symbol duration
F_sample=1/T_sample; 

num_symbols = 10;
bits = randi([0, 1],1, num_symbols); %Our data to be transmitted, 1's and 0's
x = [];
t=(0:sps*num_symbols-1)*T_sample;
for bit=bits
    pulse = zeros(1,sps);
    pulse(1) = bit*2-1; %set the first value to either a 1 or -1
    x = [x, pulse]; %add the 8 samples to the signal
end
figure(1)
stem(t, x, '.-')
grid on 


% Create our root-raised-cosine filter
filtlen = 10;      % Filter length in symbols
rolloff = 0.25;    % Filter rolloff factor
h = rcosdesign(rolloff,filtlen,sps);
num_taps=length(h);
figure(2)
plot(h, '.')
grid on

% Filter our signal, in order to apply the pulse shaping
x_shaped = conv(x, h);
t=(0:length(x_shaped)-1)*T_sample;
figure(3)
plot(t,x_shaped, 'b.-'); hold on

%% Receiver processing
x_received = conv(x_shaped, h);
t=(0:length(x_received)-1)*T_sample;
plot(t,x_received, 'r.-'); hold on

% Sample the signal at Ts, consider the delay in filter
nn=(0:num_symbols-1)*sps+(num_taps-1)+1;
X_sam=x_received(nn)
stem((nn-1)*T_sample,X_sam,'go',"filled");
% plt.grid(True)
% plt.show()

figure(4)
[pxx,f] = pspectrum(x_shaped, F_sample);
semilogy(f,pxx);grid on
