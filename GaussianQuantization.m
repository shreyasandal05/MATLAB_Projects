clear

N_obs=1e6; % No of observations
x=randn(1,1e6);
[f,xi] = ksdensity(x);
subplot(2,1,1)
plot(xi,f);hold


T1 = numerictype(1,4,3);
x1 = double(fi(x,T1));

x_vals = unique(x1);
pmf=zeros(size(x_vals));
CDF=zeros(size(x_vals));
% Finding PMF
for ii=1:length(x_vals)
    pmf(ii)=sum(x1==x_vals(ii))/N_obs;
end

% subplot(2,1,2)
stem(x_vals,pmf)
xlim([-6,6])
% Finding CDF
for ii=1:length(x_vals)
    CDF(ii)=sum(pmf(1:ii));
end
subplot(2,1,2)

stem(x_vals,CDF)
xlim([-6,6])