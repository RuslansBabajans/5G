% Practical work 1
clear all; clc; close all

Fs=100e6; % baseband clock
Fd=4000e6; % analogue sampling freq
f0=140e6; % carrier frequency
N=200; % number of symbols to transmit

% Generate user data
usrDat=kron(randi(2,1,N)*2-3,ones(1,4));

% Resampling
usrDatRsm=resample(usrDat,40,1);

% FM Modulation
t=(1:length(usrDatRsm))/Fd;
sFM=cos(2*pi*f0*t+1e-2*cumsum(usrDatRsm));

% Frequency discrimination
fr_discr=fir1(36,150e6/2e9, 'high'); % detuned contour
sFMdiscr=filter(fr_discr,1,sFM);

% AM demodulation
sFMdem=abs(sFMdiscr);

% Filter  the signal
LPF=fir1(50,f0/(Fd/2));
sFMflt=filter(LPF,1,sFMdem);
sFMrest=sFMflt-mean(sFMflt);
sFMrest=sFMrest*sqrt(mean(usrDatRsm.^2)/mean(sFMrest.^2));
%fvtool(LPF)


% Calculate FM spectra
[spectr, fr]       =win_fft(usrDatRsm, 4e9,1e4,7e3);
[spectrFM, fr]     =win_fft(sFM,       4e9,1e4,7e3);
[spectrFMdiscr, fr]=win_fft(sFMdiscr,  4e9,1e4,7e3);
[spectrFMdem, fr]  =win_fft(sFMdem,    4e9,1e4,7e3);
[spectrFMflt fr]   =win_fft(sFMrest,   4e9,1e4,7e3);


figure(1); hold on
plot(fr, 20*log10(spectr),'b.-')
plot(fr, 20*log10(spectrFM),'r.-')
plot(fr, 20*log10(spectrFMdiscr),'g.-')
plot(fr, 20*log10(spectrFMdem),'m.-')
plot(fr, 20*log10(spectrFMflt),'c.-')




figure(2); hold on
plot(usrDatRsm,'b.-')
plot(sFM,'r.-')
plot(sFMdiscr,'g.-')
plot(sFMrest,'m.-')
xlabel('Поугаи')
ylabel('Удавы')



