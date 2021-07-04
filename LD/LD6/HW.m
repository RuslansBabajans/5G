clear variables, close all, format compact

%% HW task 2

N=100000;

x=-N/2:1:N/2-1;
y=-N/2:1:N/2-1;

stan_div=1.25;
mean_val=0;

realGaus=(1/(stan_div*sqrt(2*pi)))*exp(-0.5*((x-mean_val)/(stan_div)).^2);
imagGaus=(1/(stan_div*sqrt(2*pi)))*exp(-0.5*((y-mean_val)/(stan_div)).^2);

r=sqrt(x.^2+y.^2);

ampDens_1=abs(realGaus+j*imagGaus);
phsDens_1=phase(realGaus+j*imagGaus);

ampDens_2=(r./(stan_div).^2).*exp(-(r.^2)./(2*stan_div.^2));


phsDens_2=1/(2*pi)*ones(1,N);

% figure(1)
% ksdensity(realGaus)
% hold on
% ksdensity(imagGaus)

figure(2)
ksdensity(ampDens_1)
hold on
ksdensity(ampDens_2)

figure(3)
ksdensity(phsDens_1)
hold on
ksdensity(phsDens_2)




