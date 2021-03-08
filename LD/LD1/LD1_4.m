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
%% Direct conversion receiver task

f1=170e6; % Carrier frequency 1
f2=370e6; % Carrier frequency 2

f_LO_1=f1;
f_LO_2=f2;
% Modulation

sAM_1=(2+usrDatRsm_1).*cos(2*pi*f1*t);
sAM_2=(2+usrDatRsm_2).*cos(2*pi*f2*t);

sAM=sAM_1+sAM_2;
%==============================================%
% Channel attenuation

dB_att=-4;
sAM_RX=sAM*10^(dB_att/20); %  attenuated by 4 dB
%==============================================%
% Ceparate signals using filter after antenna. 170 MHz and 370 MHz are
% image frequencies, and both will be shifted to 100 MHz intermediate
% frequency.

BP_1=fir1(50,[f1-20e6,f1+20e6]/(4e9/2),'bandpass'); % filter to separate 170 MHz carrier frequency 1
sAM_RX_1=filter(BP_1,1,sAM_RX);

BP_2=fir1(50,[f2-20e6,f2+20e6]/(4e9/2),'bandpass'); % filter to separate 370 MHz carrier frequency 2
sAM_RX_2=filter(BP_2,1,sAM_RX);
%==============================================%
% Adjusts the level of the received signal to twice the level of modulating
% signal

P_usrDatRsm_1=rms(usrDatRsm_1)^2; % s1(t) average power
P_sAM_RX_1=rms(sAM_RX_1)^2; % 170 MHz signal average power

tuned_receiver_gain_1=1; % tunable gain of the radio frequency amplifier
sAM_RX_tuned_1=sAM_RX_1; % Placeholder
Level_difference_1=10; % Placeholder

P_usrDatRsm_2=rms(usrDatRsm_2)^2; % s2(t) average power
P_sAM_RX_2=rms(sAM_RX_2)^2; % 370 MHz signal average power

tuned_receiver_gain_2=1; % tunable gain of the radio frequency amplifier
sAM_RX_tuned_2=sAM_RX_2; % Placeholder
Level_difference_2=10; % Placeholder

% Received signal amplification to adjust the power level of AM_RX to twice
% the powe level of

while Level_difference_1 > 0
    sAM_RX_tuned_1=sAM_RX_tuned_1*tuned_receiver_gain_1; % Amplify received signal
    
    P_sAM_RX_tuned_1=rms(sAM_RX_tuned_1)^2; % tuned signal average power
    
    Level_difference_1=P_usrDatRsm_1*2-P_sAM_RX_tuned_1; % tuned AM_RX power comparison to twice the level of s1(t)
    
    tuned_receiver_gain_1=tuned_receiver_gain_1+0.00001; % increase gain of the radio frequency amplifier
end

while Level_difference_2 > 0
    sAM_RX_tuned_2=sAM_RX_tuned_2*tuned_receiver_gain_2; % Amplify received signal
    
    P_sAM_RX_tuned_2=rms(sAM_RX_tuned_2)^2; % tuned signal average power
    
    Level_difference_2=P_usrDatRsm_2*2-P_sAM_RX_tuned_2; % tuned AM_RX power comparison to twice the level of s1(t)
    
    tuned_receiver_gain_2=tuned_receiver_gain_2+0.00001; % increase gain of the radio frequency amplifier
end
%==============================================%
% Mixing

LO_1=cos(2*pi*f_LO_1*t);
LO_2=cos(2*pi*f_LO_2*t-(114*pi/180)); % Addresss the phase shift after filter

sAM_1_directly_converted=sAM_RX_tuned_1.*LO_1;
sAM_2_directly_converted=sAM_RX_tuned_2.*LO_2;

% Baseband filter

conversion_filter=fir1(50,100e6/(Fd/2)); % Filter to remove the upconverted component

sAM_1_directly_converted=filter(conversion_filter,1,sAM_1_directly_converted);
sAM_1_directly_converted=sAM_1_directly_converted-mean(sAM_1_directly_converted);
sAM_1_directly_converted=sAM_1_directly_converted*sqrt(mean(usrDatRsm_1.^2)/mean(sAM_1_directly_converted.^2));

sAM_2_directly_converted=filter(conversion_filter,1,sAM_2_directly_converted);
sAM_2_directly_converted=sAM_2_directly_converted-mean(sAM_2_directly_converted);
sAM_2_directly_converted=sAM_2_directly_converted*sqrt(mean(usrDatRsm_2.^2)/mean(sAM_2_directly_converted.^2));
%==============================================%
% Spectra

[spectr_usrDatRsm_1, fr]=win_fft(usrDatRsm_1, 4e9,10^4,10^3);
[spectr_usrDatRsm_2, fr]=win_fft(usrDatRsm_2, 4e9,10^4,10^3);

[spectr_sAM_1, fr]=win_fft(sAM_1, 4e9,10^4,10^3);
[spectr_sAM_2, fr]=win_fft(sAM_2, 4e9,10^4,10^3);
[spectr_sAM, fr]=win_fft(sAM, 4e9,10^4,10^3);

[spectr_sAM_RX, fr]=win_fft(sAM_RX, 4e9,10^4,10^3);

[spectr_sAM_RX_1, fr]=win_fft(sAM_RX_1, 4e9,10^4,10^3);
[spectr_sAM_RX_2, fr]=win_fft(sAM_RX_2, 4e9,10^4,10^3);

[spectr_sAM_1_directly_converted, fr]=win_fft(sAM_1_directly_converted, 4e9,10^4,10^3);
[spectr_sAM_2_directly_converted, fr]=win_fft(sAM_2_directly_converted, 4e9,10^4,10^3);
%==============================================%
%% Plots

figure(32)
hold on
plot(t*1e6,sAM_1,'Linewidth',2)
plot(t*1e6,sAM_RX_tuned_1,'Linewidth',2)
plot(t*1e6,LO_1,'Linewidth',2)
% plot(t*1e6,sAM_1_directly_converted,'Linewidth',2)
xlim([0, 0.1])
set(gca,'fontsize',20)

figure(33)
hold on
plot(t*1e6,sAM_2,'Linewidth',2)
plot(t*1e6,sAM_RX_tuned_2,'Linewidth',2)
plot(t*1e6,LO_2,'Linewidth',2)
% plot(t*1e6,sAM_2_directly_converted,'Linewidth',2)
xlim([0, 0.1])
set(gca,'fontsize',20)

figure(1)
hold on
plot(fr*1e-6,20*log10(spectr_usrDatRsm_1))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_2))
xlim([0, 200])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("s(t)_1","s(t)_2",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(2)
subplot(4,1,1)
hold on
plot(fr*1e-6,20*log10(spectr_sAM_1))
plot(fr*1e-6,20*log10(spectr_sAM_2))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_1}","AM_{s(t)_2}",'location','northeastOutside')
grid on, grid minor
set(gca,'fontsize',20)

subplot(4,1,2)
plot(fr*1e-6,20*log10(spectr_sAM),'color','#EDB120')
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_1}+AM_{s(t)_2}",'location','northeastOutside')
grid on, grid minor
set(gca,'fontsize',20)

subplot(4,1,3)
plot(fr*1e-6,20*log10(spectr_sAM_RX_1))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_1}RX",'location','northeastOutside')
grid on, grid minor
set(gca,'fontsize',20)

subplot(4,1,4)
plot(fr*1e-6,20*log10(spectr_sAM_RX_2),'color','#D95319')
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_2}RX",'location','northeastOutside')
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(3)
subplot(2,1,1)
hold on
plot(fr*1e-6,20*log10(spectr_sAM_RX_1))
plot(fr*1e-6,20*log10(spectr_sAM_1_directly_converted))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_1} RX","AM_{s(t)_1} B.B.",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)

subplot(2,1,2)
hold on
plot(fr*1e-6,20*log10(spectr_sAM_RX_2))
plot(fr*1e-6,20*log10(spectr_sAM_2_directly_converted))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_2} RX","AM_{s(t)_2} B.B.",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)
%==============================================%

figure(4)
hold on
plot(fr*1e-6,20*log10(spectr_sAM_RX_1))
plot(fr*1e-6,20*log10(spectr_sAM_1_directly_converted))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_1))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_1} RX","AM_{s(t)_1} B.B.","s(t)_1",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)

figure(5)
hold on
plot(fr*1e-6,20*log10(spectr_sAM_RX_2))
plot(fr*1e-6,20*log10(spectr_sAM_2_directly_converted))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_2))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("AM_{s(t)_2} RX","AM_{s(t)_2} B.B.","s(t)_2",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)

figure(6)
hold on
plot(t*1e6,usrDatRsm_1)
plot(t*1e6,sAM_1_directly_converted)
ylim([-2, 2])
xlabel('t, ns')
ylabel('s(t), V')
legend("s(t)_1","s(t)_1 Demodulated",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)

figure(7)
hold on
plot(t*1e6,usrDatRsm_2)
plot(t*1e6,sAM_2_directly_converted)
ylim([-2, 2])
xlabel('t, ns')
ylabel('s(t), V')
legend("s(t)_2","s(t)_2 Demodulated",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)



