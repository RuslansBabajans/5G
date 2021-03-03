%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters

Fs=100e6; % 100 MHz baseband clock
Fd=4000e6; % analogue sampling freq.
N=200; % number of symbols to transmit


% Generate user data
usrDat_1=kron(randi(2,1,N)*2-3,ones(1,4)); % ones(1,4) is 4 samples per symbol,  Fs/4= 25 Mbaud
usrDat_2=kron(randi(2,1,N)*2-3,ones(1,4));

% Resampling
usrDatRsm_1=resample(usrDat_1,40,1);
usrDatRsm_2=resample(usrDat_2,40,1);

% Time vector
t=(1:length(usrDatRsm_1))/Fd;
%==============================================%
%% Superheterodyne receiver task

f1=170e6; % Carrier frequency 1
f2=370e6; % Carrier frequency 2

f_if=100e6; % Intemideate frequency
f_LO=270e6; % f_LO_1 = f1 + f_i  and  f_LO_2 = f2 - f_if are the same. Mixing with f_LO will invert spectrum of f1

% Modulation

sFM_1=cos(2*pi*f1*t+1e-2*cumsum(usrDatRsm_1));
sFM_2=cos(2*pi*f2*t+1e-2*cumsum(usrDatRsm_2));

sFM=sFM_1+sFM_2;
%==============================================%
% Channel attenuation

dB_att=-4;
sFM_RX=sFM*10^(dB_att/20); %  attenuated by 4 dB
%==============================================%
% Adjusts the level of the received signal to twice the level of modulating
% signal

P_usrDatRsm_1=rms(usrDatRsm_1)^2; % s1(t) average power
P_sFM_RX=rms(sFM_RX)^2; % AM_RX average power

tuned_receiver_gain=1; % tunable gain of the radio frequency amplifier

sFM_RX_tuned=sFM_RX; % Placeholder
Level_difference=10; % Placeholder 

% Received signal amplification to adjust the power level of AM_RX to twice
% the powe level of s(t)

while Level_difference > 0
sFM_RX_tuned=sFM_RX_tuned*tuned_receiver_gain; % Amplify received signal

P_sFM_RX_tuned=rms(sFM_RX_tuned)^2; % tuned signal average power

Level_difference=P_usrDatRsm_1*2-P_sFM_RX_tuned; % tuned AM_RX power comparison to twice the level of s1(t) 

tuned_receiver_gain=tuned_receiver_gain+0.00000001; % increase gain of the radio frequency amplifier
end
%==============================================%
% Mixing

LO=1*cos(2*pi*f_LO*t);

sFM_at_IF=sFM_RX_tuned.*LO;

%==============================================%
% Demodulation using diode

sFMdem=abs(sFM_at_IF);
%==============================================%
% Add white Gaussian noise to the signal produced by the diode

standard_deviation=0.1; % Standard deviation
Gaussian_noise=randn(size(sFMdem));

std_check=std(standard_deviation*Gaussian_noise); % Check Standard deviation, std_check=0.09981

sAMdem_and_Gaussian_noise=sFMdem+standard_deviation*Gaussian_noise;
%==============================================%
% Filter  the signal

% LPF=fir1(50,f0/(Fd/2));
% sAMflt=filter(LPF,1,sAMdem_and_Gaussian_noise);
%==============================================%
% Spectra 
% [spectr, fr]=win_fft(usrDatRsm, 4e9,10^4,10^3);
% [spectrAM_TX, fr]=win_fft(sAM_TX, 4e9,10^4,10^3);
% [spectrAM_RX, fr]=win_fft(sAM_RX, 4e9,10^4,10^3);
% [spectrAMdem, fr]=win_fft(sAMdem, 4e9,10^4,10^3);
% [spectrsAMdem_and_Gaussian_noise, fr]=win_fft(sAMdem_and_Gaussian_noise, 4e9,10^4,10^3);
% [spectrAMflt, fr]=win_fft(sAMflt, 4e9,10^4,10^3);
%==============================================%
%% Plots
% 
% figure(1)
% hold on
% plot(t*1e6,sAM_TX)
% plot(t*1e6,sAM_RX)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("AM_{TX}","AM_{RX}")
% grid on, grid minor
% set(gca,'fontsize',12)
% %==============================================%
% 
% figure(2)
% hold on
% plot(t*1e6,sAM_RX_tuned)
% plot(t*1e6,sAM_RX)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("AM_{RX} tuned","AM_{RX}")
% grid on, grid minor
% set(gca,'fontsize',12)
% %==============================================%
% 
% figure(3)
% 
% subplot(4,1,1)
% plot(t,usrDatRsm)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("usrDat",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% 
% subplot(4,1,2)
% plot(t*1e6,sAMdem)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("AMdem",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% 
% subplot(4,1,3)
% plot(t*1e6,sAMdem_and_Gaussian_noise)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("AMdem+Gaussian noise",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% 
% subplot(4,1,4)
% plot(t*1e6,sAMflt)
% 
% xlabel('t, ns')
% ylabel('s(t), V')
% legend("AMflt",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% %==============================================%
% 
% figure(4)
% plot(fr,20*log10(spectrAM_TX))
% hold on
% plot(fr,20*log10(spectrAM_RX)-4)
% 
% xlabel('f, Hz')
% ylabel('s(f), dB')
% legend("AM_{TX}","AM_{RX}",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% %==============================================%
% 
% figure(5)
% plot(fr, 20*log10(spectr))
% hold on
% plot(fr,20*log10(spectrAMdem))
% plot(fr,20*log10(spectrsAMdem_and_Gaussian_noise))
% plot(fr,20*log10(spectrAMflt))
% 
% xlabel('f, Hz')
% ylabel('s(f), dB')
% legend("usrDat","AMdem","AMdem+Gaussian noise","AMflt",'location','northeast')
% grid on, grid minor
% set(gca,'fontsize',12)
% %==============================================%




