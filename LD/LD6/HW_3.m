clear variables, close all, format compact

%% HW task 2

N=100000;

x=linspace(0,3,N);
y=linspace(0,3,N);

stan_div=1.25;
mean_val_1=2;
mean_val_2=3;

realGaus=stan_div*randn(N,1)+mean_val_1;
imagGaus=stan_div*randn(N,1)+mean_val_2;

r=sqrt(x.^2+y.^2);

ampDens_1=abs(realGaus+j*imagGaus);


% z=((mean_val_1*mean_val_2)/(2*stan_div^2)).*r;
% ampDens_2=(r./(stan_div^2)).*exp(-1.*(r.^2+mean_val_1^2+mean_val_2^2)./(2*stan_div^2)).*besselj(0,z);


figure(1)
histogram(realGaus)
hold on
histogram(imagGaus)

figure(2)
histogram(ampDens_1)
% hold on
% plot(r,ampDens_2*1e5)








