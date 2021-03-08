% Practical work 2
clear all; clc; close all

Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
N=1000; % number of symbols to transmit

ALPH=-3:2:3;

% Generate user data
usrDatI=kron(ALPH(randi(4,1,N)),[1 0 0 0]);
usrDatQ=kron(ALPH(randi(4,1,N)),[1 0 0 0]);

% Pulse-sahping filter
%firTx=firrcos(46,0.25,0.25,2,'rolloff','normal');
firTx=firrcos(46,0.25,0.25,2,'rolloff', 'sqrt');
%fvtool(firTx,1,firTx2,1,conv(firTx2,firTx2),1)
% return

usrDatFltI=filter(firTx, 1, usrDatI); usrDatFltI=usrDatFltI(46:end);
usrDatFltQ=filter(firTx, 1, usrDatQ); usrDatFltQ=usrDatFltQ(46:end);

% figure(1); hold on
% plot(usrDatI(1:100), 'b.-')
% plot(usrDatFltI(1:100), 'r.-')
% return


% Resampling
usrDatRsmI=resample(usrDatFltI,40,1);
usrDatRsmQ=resample(usrDatFltQ,40,1);

% Modulation - I*cos-Q*sin
t=(1:length(usrDatRsmI))/Fd;
sAM=usrDatRsmI.*cos(2*pi*f0*t)-usrDatRsmQ.*sin(2*pi*f0*t);


% QAM demodulation
sQAMdemI= sAM.*cos(2*pi*f0*t);
sQAMdemQ=-sAM.*sin(2*pi*f0*t);

% Filter  the signal
LPF=fir1(100,f0/(Fd/2));
sAMfltI=2*filter(LPF,1,sQAMdemI); sAMfltI=sAMfltI(161:end);
sAMfltQ=2*filter(LPF,1,sQAMdemQ); sAMfltQ=sAMfltQ(161:end);
%fvtool(LPF)

% Resampling
usrDatRsm2I=resample(sAMfltI(11:end),1,40);
usrDatRsm2Q=resample(sAMfltQ(11:end),1,40);

% figure(1); hold on
% plot(usrDatFltQ,'b.-')
% plot(usrDatRsm2Q,'r.-')

usrDatRsm2I=filter(firTx, 1, usrDatRsm2I); usrDatRsm2I=usrDatRsm2I(60:end);
usrDatRsm2Q=filter(firTx, 1, usrDatRsm2Q); usrDatRsm2Q=usrDatRsm2Q(60:end);


% eyediagram(usrDatRsm2Q(4:end),16)
% return


plot(usrDatRsm2I(4:4:end),usrDatRsm2Q(4:4:end),'b.')
return

% Calculate AM spectra
[spectr, fr]=win_fft(usrDatRsmI, 4e9,10^4,10^3);
[spectrAM, fr]=win_fft(sAM, 4e9,10^4,10^3);
[spectrQAMdem, fr]=win_fft(sQAMdemQ, 4e9,10^4,10^3);
[spectrQAMflt, fr]=win_fft(sAMfltQ, 4e9,10^4,10^3);


% figure(1); hold on; grid on
% plot(fr, 20*log10(spectr),'b.-')
% plot(fr, 20*log10(spectrAM),'r.-')
% plot(fr, 20*log10(spectrQAMdem),'g.-')
% plot(fr, 20*log10(spectrQAMflt),'m.-')
% return

% figure(1); hold on
% plot(usrDatRsm,'b.-')
% plot(sAM,'r.-')
% plot(sAMdem,'g.-')
% plot(sAMflt,'m.-')

