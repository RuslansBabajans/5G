clear variables, close all, format compact

%% HW task 2

N=100000;

stan_div=1.25;
mean_val_1=2;
mean_val_2=3;

realGaus=stan_div*randn(N,1)+mean_val_1;
imagGaus=stan_div*randn(N,1)+mean_val_2;

r=linspace(0,8,N);

ampDens_1=abs(realGaus+j*imagGaus);

mean_total=sqrt(mean_val_1^2+mean_val_2^2);
z=((mean_total*r)/(stan_div^2));
ampDens_2=(r./(stan_div^2)).*besseli(0,z).*exp(-1*(r.^2+mean_val_1^2+mean_val_2^2)./(2*stan_div^2));


figure(1)
histogram(realGaus)
hold on
histogram(imagGaus)

figure(2)
histogram(ampDens_1)
hold on
plot(r,ampDens_2*1e4,'LineWidth',2)

% figure(3)
% plot(r,ampDens_2*1e4)








