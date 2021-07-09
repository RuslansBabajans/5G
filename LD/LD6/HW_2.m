clear variables, close all, format compact

%% HW task 2

N=100000;

x=linspace(0,5,N);
y=linspace(0,5,N);

stan_div=1.25;
mean_val=0;

realGaus=stan_div*randn(N,1)+mean_val;
imagGaus=stan_div*randn(N,1)+mean_val;

r=sqrt(x.^2+y.^2);
% phi=atan(y./x);

ampDens_1=abs(realGaus+j*imagGaus);
phsDens_1=angle(realGaus+j*imagGaus);

ampDens_2=(r./((stan_div)^2)).*exp(-(r.^2)/(2*stan_div^2));
phsDens_2=1600*1/(2*pi)*ones(1,100); %rand(N,1)*2*pi-pi;

phi=linspace(-pi,pi,100);

figure(1)
histogram(realGaus)
hold on
histogram(imagGaus)

figure(2)
[d]=histogram(ampDens_1,100)
hold on
plot(r,ampDens_2*5e3,'LineWidth',2)

figure(3)
histogram(phsDens_1)
hold on
plot(phi,phsDens_2)







