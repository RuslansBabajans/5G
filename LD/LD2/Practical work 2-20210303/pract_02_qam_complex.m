% Practical work 2
clear all; clc; close all

Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
N=1000; % number of symbols to transmit

ALPH=-3:2:3; ALPH=repmat(-3:2:3,4,1)+1j*repmat(ALPH.',1,4);

% Generate user data
usrDat=kron(ALPH(randi(numel(ALPH),1,N)),[1 0 0 0]);


% Pulse-sahping filter
%firTx=firrcos(46,0.25,0.25,2,'rolloff','normal');
firTx=firrcos(46,0.25,0.25,2,'rolloff', 'sqrt');


usrDatFlt=filter(firTx, 1, usrDat); usrDatFlt=usrDatFlt(46:end);


t=(1:length(usrDatFlt))/Fs;
usrChan=usrDatFlt.*exp(j*2*pi*0e3*t);
SNR=20; % signal-to-noise ratio in dB
usrChan=usrChan+(randn(size(usrChan))+1j*randn(size(usrChan)))*sqrt(mean(abs(usrChan).^2)/2)*10^(-SNR/20);


usrDatRsm2=filter(firTx, 1, usrChan); usrDatRsm2=usrDatRsm2(60:end);



% eyediagram(usrDatRsm2Q(3:end),16)
% return

figure(1)
plot(usrDatRsm2(3:4:end),'b.')

% Calculate AM spectra
[spectr, fr]=win_fft(usrDatRsm2, 2,1000,10);
%[spectrAM, fr]=win_fft(sAM, 4e9,10^4,10^3);
%[spectrQAMdem, fr]=win_fft(sQAMdemQ, 4e9,10^4,10^3);
% [spectrQAMflt, fr]=win_fft(sAMfltQ, 4e9,10^4,10^3);


figure(2); hold on; grid on
plot(fr, 20*log10(spectr),'b.-')
% plot(fr, 20*log10(spectrAM),'r.-')
% plot(fr, 20*log10(spectrQAMdem),'g.-')
% plot(fr, 20*log10(spectrQAMflt),'m.-')


% figure(3); hold on
% plot(usrDatRsm,'b.-')
% plot(sAM,'r.-')
% plot(sAMdem,'g.-')
% plot(sAMflt,'m.-')

