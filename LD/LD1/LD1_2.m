%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters

Fs=100e6; % 100 MHz baseband clock
Fd=4000e6; % analogue sampling freq.
N=200; % number of symbols to transmit


% Generate user data usrDat=kron(randi(2,1,N)*2-3,ones(1,4));
usrDat=kron(randi(2,1,N)*2-3,ones(1,4)); % ones(1,4) is 4 samples per symbol,  Fs/4= 25 Mbaud

% Resampling
usrDatRsm=resample(usrDat,40,1);

% Time vector
t=(1:length(usrDatRsm))/Fd;
%==============================================%
%% Tuned radio frequency receiver task

f0=140e6; % Carrier frequency

% Modulation
sAM_TX=(2+usrDatRsm).*cos(2*pi*f0*t);
%==============================================%
% Channel attenuation

dB_att=-4;
sAM_RX=sAM_TX*10^(dB_att/20); %  sAM_RX= sAM_TX attenuated by 4 dB
attenuation_check=20*log10(sAM_RX./sAM_TX); % attenuation_check= -4.000000 dB
%==============================================%
% Adjusts the level of the received signal to twice the level of modulating
% signal

P_usrDatRsm=rms(usrDatRsm)^2; % s(t) average power
P_sAM_RX=rms(sAM_RX)^2; % AM_RX average power

tuned_receiver_gain=1; % tunable gain of the radio frequency amplifier

sAM_RX_tuned=sAM_RX; % Placeholder
Level_difference=10; % Placeholder

% Received signal amplification to adjust the power level of AM_RX to twice
% the powe level of s(t)

while Level_difference > 0
    sAM_RX_tuned=sAM_RX*tuned_receiver_gain; % Amplify received signal
    
    P_sAM_RX_tuned=rms(sAM_RX_tuned)^2; % tuned signal average power
    
    Level_difference=P_usrDatRsm*2-P_sAM_RX_tuned; % tuned AM_RX power comparison to twice the level of s(t)
    
    tuned_receiver_gain=tuned_receiver_gain+0.00001; % increase gain of the radio frequency amplifier
end

%==============================================%
% AM demodulation using crystal receiver (diode)

sAMdem=abs(sAM_RX_tuned);
%==============================================%
% Add white Gaussian noise to the signal produced by the diode

standard_deviation=0.1; % Standard deviation
Gaussian_noise=randn(size(sAMdem));

std_check=std(standard_deviation*Gaussian_noise); % Check Standard deviation, std_check=0.09981

sAMdem_and_Gaussian_noise=sAMdem+standard_deviation*Gaussian_noise;
%==============================================%
% Filter  the signal

LPF=fir1(50,f0/(Fd/2));
sAMflt=filter(LPF,1,sAMdem_and_Gaussian_noise);
%==============================================%
% Spectra
[spectr, fr]=win_fft(usrDatRsm, 4e9,10^4,10^3);
[spectrAM_TX, fr]=win_fft(sAM_TX, 4e9,10^4,10^3);
[spectrAM_RX, fr]=win_fft(sAM_RX, 4e9,10^4,10^3);
[spectrAMdem, fr]=win_fft(sAMdem, 4e9,10^4,10^3);
[spectrsAMdem_and_Gaussian_noise, fr]=win_fft(sAMdem_and_Gaussian_noise, 4e9,10^4,10^3);
[spectrAMflt, fr]=win_fft(sAMflt, 4e9,10^4,10^3);
%==============================================%
%% Plots

figure(1)
hold on
plot(t*1e6,sAM_TX)
plot(t*1e6,sAM_RX)

xlabel('t, ns')
ylabel('s(t), V')
legend("AM_{TX}","AM_{RX}")
grid on, grid minor
set(gca,'fontsize',20)

figure(5)
hold on
plot(t*1e6,sAM_RX_tuned)
plot(t*1e6,sAM_RX)

xlabel('t, ns')
ylabel('s(t), V')
legend("AM_{RX}","AM_{RX} tuned")
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(2)

subplot(2,1,1)
plot(t*1e6,usrDatRsm)

xlabel('t, ns')
ylabel('s(t), V')
title("Modulating signal")
grid on, grid minor
set(gca,'fontsize',20)

subplot(2,1,2)
plot(t*1e6,sAMdem)

xlabel('t, ns')
ylabel('s(t), V')
title("Diode output")
grid on, grid minor
set(gca,'fontsize',20)


figure(8)

subplot(2,1,1)
plot(t*1e6,sAMdem_and_Gaussian_noise)

xlabel('t, ns')
ylabel('s(t), V')
title("Diode output+Gaussian noise")
grid on, grid minor
set(gca,'fontsize',20)

subplot(2,1,2)
plot(t*1e6,sAMflt)

xlabel('t, ns')
ylabel('s(t), V')
title("Demodulated signal")
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(3)
plot(fr*1e-6,20*log10(spectrAM_TX))
hold on
plot(fr*1e-6,20*log10(spectrAM_RX)-4)
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{TX}","AM_{RX}",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(4)
plot(fr*1e-6, 20*log10(spectr))
hold on
plot(fr*1e-6,20*log10(spectrAMdem))
plot(fr*1e-6,20*log10(spectrsAMdem_and_Gaussian_noise))
plot(fr*1e-6,20*log10(spectrAMflt))

xlabel('f, MHz')
ylabel('s(f), dB')
legend("Modulating signal","Diode output","Diode output+Gaussian noise","Demodulated signal",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%
