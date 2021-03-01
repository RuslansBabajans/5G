% Practical work N1

%% Clear variables and close figures
format long
clear variables
close all
%=============================%

%% Parameters

Fs=100e6; %100Mhz baseband clock
Fd=4e9; % 4GHz analogue sampling frequency
f0=140e6; % carrier frequency
N=200; % number of symbols
usrDat=kron(randi([-1,1],1,N)-1,ones(1,4));

% Resampling
usrDatrsm=resample(usrDat,40,1);

% Modulation
t=(1:length(usrDatrsm))/Fd;
sAM=(2+usrDatrsm).*cos(2*pi*f0*t);

% figure(1)
% plot(usrDat,'b.-')

figure(2)
plot(usrDatrsm,'b.-')
hold on
plot(sAM,'r.-')

% figure(3)
% plot(sAM,'r.-')

%=============================%

%% Spectrum

spectr_AM=fft(sAM);
spectr_BB=fft(usrDatrsm);
amp_spectr_AM_db=20*log10(abs(spectr_AM));
amp_spectr_BB_db=20*log10(abs(spectr_BB));

figure(4)
plot(amp_spectr_AM_db,'r.-')
hold on
plot(amp_spectr_BB_db,'b.-')

%=============================%
%% AM demodulation

asMdem=abs(sAM);

amp_spectr_asMdem_db=20*log10(abs(fft(asMdem)));


figure(5)
plot(amp_spectr_AM_db,'r.-')
hold on
plot(amp_spectr_BB_db,'b.-')
plot(amp_spectr_asMdem_db,'g.-')

%=============================%
%% Filter

LPF=fir1(50,f0/(Fd/2));

sAM_flt=filter(LPF,1,asMdem);
amp_sAM_flt_db=20*log10(abs(fft(sAM_flt)));

% fvtool(LPF)
% return

figure(10)
plot(amp_spectr_AM_db,'r.-')
hold on
plot(amp_spectr_BB_db,'b.-')
plot(amp_spectr_asMdem_db,'g.-')
plot(amp_sAM_flt_db,'m.-')

figure(11)
plot(usrDatrsm,'r.-')
hold on
plot(sAM_flt,'b.-')




