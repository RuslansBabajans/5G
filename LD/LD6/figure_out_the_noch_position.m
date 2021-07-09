close all; clear variables; format short

f0=1e9;
tau=5e-9;

% w=linspace(0,2*pi,10000);
% f=w/(2*pi);
f=linspace(0,1e6,1e6);

A=1;
B=1;
amp=abs((-A-B*exp(-j*2*pi*f0*tau))./(j*2*pi*f));

figure(1)
plot(f,amp)