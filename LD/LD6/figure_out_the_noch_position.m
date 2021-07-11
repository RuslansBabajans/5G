close all; clear variables; format short

f0=1e9;
tau=5e-9;

% w=linspace(0,2*pi,10000);
% f=w/(2*pi);
f=linspace(0,1/tau,1e6);

DP=10^(-3/20);
A=1/sqrt(1+(1-DP)^2);
B=(1-DP)/sqrt(1+(1-DP)^2);


w=2*pi*tau*(f0+f);

amp=abs(A+B*exp(-j*2*pi*tau*(f0+f)));

[minimum, minIndex]=min(amp);

figure(1)
plot(f/1e8-1,amp)
hold on
plot(f(minIndex)/1e8-1,minimum,'o')





%% Plot the dependence of the position of this minimum on the time shift 
%% τ ∈ (4.5, 5.5) ns.

tau=linspace(4.51e-9, 5.49e-9,100);

minimum=zeros(1,length(tau));
minIndex=zeros(1,length(tau));


for n=1:1:length(tau)
f=linspace(0,1/tau(n),1e6);

DP=10^(-3/20);
A=1/sqrt(1+(1-DP)^2);
B=(1-DP)/sqrt(1+(1-DP)^2);



amp=abs(A+B*exp(-j*2*pi*tau(n)*(f0+f)));

[minimum(n), minIndex(n)]=min(amp);
end

figure(2)
plot(tau,minimum)

figure(2)
plot(tau,f(minIndex)/1e8-1)
