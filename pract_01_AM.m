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

% Modulation
t=(1:length(usrDatRsm))/Fd;
sAM=(2+usrDatRsm).*cos(2*pi*f0*t);


% AM demodulation
sAMdem=abs(sAM);

% Filter  the signal
LPF=fir1(50,f0/(Fd/2));
sAMflt=filter(LPF,1,sAMdem);
%fvtool(LPF)


% Calculate AM spectra
[spectr, fr]=win_fft(usrDatRsm, 4e9,10^4,10^3);
[spectrAM, fr]=win_fft(sAM, 4e9,10^4,10^3);
[spectrAMdem, fr]=win_fft(sAMdem, 4e9,10^4,10^3);
[spectrAMflt, fr]=win_fft(sAMflt, 4e9,10^4,10^3);


figure(1); hold on; grid on
plot(fr, 20*log10(spectr),'b.-')
plot(fr, 20*log10(spectrAM),'r.-')
plot(fr, 20*log10(spectrAMdem),'g.-')
plot(fr, 20*log10(spectrAMflt),'m.-')
return

figure(2); hold on
plot(usrDatRsm,'b.-')
plot(sAM,'r.-')
plot(sAMdem,'g.-')
plot(sAMflt,'m.-')

