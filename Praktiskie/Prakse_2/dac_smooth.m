clear all; clc

% DAC smoothing filter demonstration


ALPH=-3:2:3;
N=1000;

fir=firrcos(36,0.25,0.25,2,'rolloff');
%fvtool(fir)

fir_sm=fir1(50,0.25);

x=kron(ALPH(randi(4,1,N)),[1 0 0 0]);
y1=filter(fir,1,x);


y2=kron(y1,[1 0 0]);

y3=filter([1 1 1],1,y2);


y4=filter(fir_sm,1,y3);

[y1sp f1]=win_fft(y1-1j*1e-7,2,1000,10);
[y2sp f2]=win_fft(y2-1j*1e-7,6,1000,10);
[y3sp f3]=win_fft(y3-1j*1e-7,6,1000,10);
[y4sp f4]=win_fft(y4-1j*1e-7,6,1000,10);

figure(1); hold on
subplot(3,2,2); plot(f1, 20*log10(y1sp),'b.-')
subplot(3,2,4); plot(f2, 20*log10(y2sp),'b.-'); hold on; plot(f3, 20*log10(y3sp),'r.-')
subplot(3,2,6); plot(f4, 20*log10(y4sp),'b.-')
subplot(3,2,1); plot(y1(100:140),'b.-'), xlim([0,40])
subplot(3,2,3); plot(y2(300:420),'b.-'); hold on; plot(y3(300:420),'r.-'), xlim([0,120])
subplot(3,2,5); plot(y4(300+25:420+25),'b.-'), xlim([0,120])



