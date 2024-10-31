clear 

Ts=1e-3;% 1 ms
sps = 8;
T_sample=Ts/sps;          % sps samples in each symbol duration
F_sample=1/T_sample; 
map=[-3,-1,1,3];

num_symbols = 20;
bits = randi([0, 1],1, num_symbols*2); %Our data to be transmitted, 1's and 0's
x = [];
t=(0:sps*num_symbols-1)*T_sample;
for jj=1:num_symbols
    pulse = zeros(1,sps);
    bb = bits((jj-1)*2+1:(jj-1)*2+2) 
    pulse(1)=map(bit2int(bb',2)+1);

    x = [x, pulse]; %add the 8 samples to the signal
end
figure(1)
stem(t, x, '.-')
grid on 


% Create our raised-cosine filter
num_taps = 101; % approximately 101/8 symbols duration
beta = 0.35;
n = (0:num_taps-1)-(num_taps-1)/2;
t=n*T_sample;
h = sinc(t/Ts) .* cos(pi*beta*t/Ts) ./ (1 - (2*beta*t/Ts).^2);
% h=ones(1,length(n));
figure(2)
plot(t, h, '.')
grid on

% Filter our signal, in order to apply the pulse shaping
x_shaped = conv(x, h);
t=(0:length(x_shaped)-1)*T_sample;
figure(3)
plot(t,x_shaped, '.-'); hold on

% Sample the signal at Ts, consider the delay in filter
nn=(0:num_symbols-1)*sps+(num_taps-1)/2+1;
X_sam=x_shaped(nn)
stem((nn-1)*T_sample,X_sam);
% plt.grid(True)
% plt.show()

figure(4)
[pxx,f] = pspectrum(x_shaped, F_sample);
semilogy(f,pxx);grid on
