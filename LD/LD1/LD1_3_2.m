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
% Ceparate signals using filter after antenna. 170 MHz and 370 MHz are
% image frequencies, and both will be shifted to 100 MHz intermediate frequency. 


BP_1=fir1(50,[f1-40e6,f1+40e6]/(4e9/2),'bandpass'); % filter to separate 170 MHz carrier frequency 1
sFM_RX_1=filter(BP_1,1,sFM_RX);


BP_2=fir1(50,[f2-40e6,f2+40e6]/(4e9/2),'bandpass'); % filter to separate 370 MHz carrier frequency 2
sFM_RX_2=filter(BP_2,1,sFM_RX);
%==============================================%
% Adjusts the level of the received signal to twice the level of modulating
% signal

P_usrDatRsm_1=rms(usrDatRsm_1)^2; % s1(t) average power
P_sFM_RX_1=rms(sFM_RX_1)^2; % 170 MHz signal average power

tuned_receiver_gain_1=1; % tunable gain of the radio frequency amplifier
sFM_RX_tuned_1=sFM_RX_1; % Placeholder
Level_difference_1=10; % Placeholder

P_usrDatRsm_2=rms(usrDatRsm_2)^2; % s2(t) average power
P_sFM_RX_2=rms(sFM_RX_2)^2; % 370 MHz signal average power

tuned_receiver_gain_2=1; % tunable gain of the radio frequency amplifier
sFM_RX_tuned_2=sFM_RX_2; % Placeholder
Level_difference_2=10; % Placeholder

% Received signal amplification to adjust the power level of AM_RX to twice
% the powe level of s(t)

while Level_difference_1 > 0
sFM_RX_tuned_1=sFM_RX_tuned_1*tuned_receiver_gain_1; % Amplify received signal

P_sFM_RX_tuned_1=rms(sFM_RX_tuned_1)^2; % tuned signal average power

Level_difference_1=P_usrDatRsm_1*2-P_sFM_RX_tuned_1; % tuned AM_RX power comparison to twice the level of s1(t) 

tuned_receiver_gain_1=tuned_receiver_gain_1+0.00001; % increase gain of the radio frequency amplifier
end

while Level_difference_2 > 0
sFM_RX_tuned_2=sFM_RX_tuned_2*tuned_receiver_gain_2; % Amplify received signal

P_sFM_RX_tuned_2=rms(sFM_RX_tuned_2)^2; % tuned signal average power

Level_difference_2=P_usrDatRsm_2*2-P_sFM_RX_tuned_2; % tuned AM_RX power comparison to twice the level of s1(t) 

tuned_receiver_gain_2=tuned_receiver_gain_2+0.00001; % increase gain of the radio frequency amplifier
end
%==============================================%
% Mixing

LO=cos(2*pi*f_LO*t);

sFM_1_at_IF=sFM_RX_tuned_1.*LO;
sFM_2_at_IF=sFM_RX_tuned_2.*LO;

% Intermediate frequency filter

IF_filter=fir1(50,f_if/(Fd/4)); % If filter to remove the upconverted component

sFM_1_at_IF=filter(IF_filter,1,sFM_1_at_IF);
sFM_2_at_IF=filter(IF_filter,1,sFM_2_at_IF);

%==============================================%
% Demodulation using diode

fr_discr=fir1(50,(110e6)/(Fd/2), 'high');
sFM_1_at_IF=filter(fr_discr,1,sFM_1_at_IF);
sFM_2_at_IF=filter(fr_discr,1,sFM_2_at_IF);

sFMdem_1=abs(sFM_1_at_IF);
sFMdem_2=abs(sFM_2_at_IF);

LPF=fir1(50,f_if/(Fd)); % Filter after the diode

sFMdem_1=filter(LPF,1,sFMdem_1);
sFMdem_1=sFMdem_1-mean(sFMdem_1);
sFMdem_1=sFMdem_1*sqrt(mean(usrDatRsm_1.^2)/mean(sFMdem_1.^2));


sFMdem_2=filter(LPF,1,sFMdem_2);
sFMdem_2=sFMdem_2-mean(sFMdem_2);
sFMdem_2=sFMdem_2*sqrt(mean(usrDatRsm_2.^2)/mean(sFMdem_2.^2));
%==============================================%
% Spectra 

[spectr_usrDatRsm_1, fr]=win_fft(usrDatRsm_1, 4e9,10^4,10^3);
[spectr_usrDatRsm_2, fr]=win_fft(usrDatRsm_2, 4e9,10^4,10^3);

[spectr_sFM_1, fr]=win_fft(sFM_1, 4e9,10^4,10^3);
[spectr_sFM_2, fr]=win_fft(sFM_2, 4e9,10^4,10^3);
[spectr_sFM, fr]=win_fft(sFM, 4e9,10^4,10^3);

[spectr_sFM_RX, fr]=win_fft(sFM_RX, 4e9,10^4,10^3);

[spectr_sFM_RX_1, fr]=win_fft(sFM_RX_1, 4e9,10^4,10^3);
[spectr_sFM_RX_2, fr]=win_fft(sFM_RX_2, 4e9,10^4,10^3);

[spectr_sFM_1_at_IF, fr]=win_fft(sFM_1_at_IF, 4e9,10^4,10^3);
[spectr_sFM_2_at_IF, fr]=win_fft(sFM_2_at_IF, 4e9,10^4,10^3);

[spectr_sFMdem_1, fr]=win_fft(sFMdem_1, 4e9,10^4,10^3);
[spectr_sFMdem_2, fr]=win_fft(sFMdem_2, 4e9,10^4,10^3);

%==============================================%
%% Plots

figure(1)
hold on
plot(fr*1e-6,20*log10(spectr_usrDatRsm_1))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_2))
xlim([0, 200])
%==============================================%

figure(2)
subplot(4,1,1)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_1))
plot(fr*1e-6,20*log10(spectr_sFM_2))
xlim([0, 500])

subplot(4,1,2)
plot(fr*1e-6,20*log10(spectr_sFM),'color','#EDB120')
xlim([0, 500])

subplot(4,1,3)
plot(fr*1e-6,20*log10(spectr_sFM_RX_1))
xlim([0, 500])

subplot(4,1,4)
plot(fr*1e-6,20*log10(spectr_sFM_RX_2),'color','#D95319')
xlim([0, 500])
%==============================================%

figure(3)
subplot(2,1,1)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_RX_1))
plot(fr*1e-6,20*log10(spectr_sFM_1_at_IF)) % Note that the spectrum is flipped
xlim([0, 500])

subplot(2,1,2)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_RX_2))
plot(fr*1e-6,20*log10(spectr_sFM_2_at_IF))
xlim([0, 500])
%==============================================%

figure(4)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_RX_1))
plot(fr*1e-6,20*log10(spectr_sFM_1_at_IF))
plot(fr*1e-6,20*log10(spectr_sFMdem_1))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_1))
xlim([0, 500])

figure(5)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_RX_2))
plot(fr*1e-6,20*log10(spectr_sFM_2_at_IF))
plot(fr*1e-6,20*log10(spectr_sFMdem_2))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_2))
xlim([0, 500])

figure(10)
hold on
plot(t*1e6,usrDatRsm_1)
plot(t*1e6,sFMdem_1)

figure(11)
hold on
plot(t*1e6,usrDatRsm_2)
plot(t*1e6,sFMdem_2)
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

% figure(4)
% plot(fr,20*log10(spectrAM_TX))
% hold on
% plot(fr,20*log10(spectrAM_RX)-4)

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




