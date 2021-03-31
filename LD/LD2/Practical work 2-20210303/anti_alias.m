clear all; clc

% ADC anti-aliasing filter demonstration


ALPH=-3:2:3;
N=1000;

fir=firrcos(66,0.125,0.25,2,'rolloff');
% fvtool(fir)

fir_sm=fir1(50,0.25);

x1=kron(ALPH(randi(4,1,N)),[1 zeros(1,7)]);
y1=filter(fir,1,x1);


fir_int=firrcos(166,0.125/2,0.25,2,'rolloff');
x2=kron(ALPH(randi(4,1,N)),[1 zeros(1,15)]);

y2=filter(fir_int,1,x2);
%y2=y2.*sin(2*pi*(0.5-0.125/4)*(1:length(y2))-pi/4);

y3=y1+2*y2(1:length(y1)).*sin(2*pi*(0.5-0.125/4)*(1:length(y1))-pi/4);
y4=y1+2*y2(1:length(y1)).*sin(2*pi*(0.5-0.2)*(1:length(y1))-pi/4);

% y2=kron(y1,[1 0 0]);
% 
% y3=filter([1 1 1],1,y2);
% 
% 
% y4=filter(fir_sm,1,y3);

[y3sp f1]=win_fft(y3,2,1000,10);
[y4sp f2]=win_fft(y4,2,1000,10);
%[y3sp f3]=win_fft(y3,6,1000,10);
%[y4sp f4]=win_fft(y4,6,1000,10);

figure(1); hold on
plot(f1, 20*log10(y3sp),'b.-')
plot(f1, 20*log10(y4sp),'r.-')


y5=y3(1:2:end);
y6=y4(1:2:end);
[y5sp f1]=win_fft(y5,2,1000,10);
[y6sp f2]=win_fft(y6,2,1000,10);


figure(2); hold on
plot(f1, 20*log10(y5sp),'b.-')
plot(f1, 20*log10(y6sp),'r.-')




