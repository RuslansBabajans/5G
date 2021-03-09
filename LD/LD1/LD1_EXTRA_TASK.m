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

% Resampling
usrDatRsm_1=resample(usrDat_1,40,1);

% Time vector
t=(1:length(usrDatRsm_1))/Fd;

%==============================================%
%% Extra task, demodulate FM using QAM demodulator
f1=170e6; % Carrier frequency 1
sFM_1=cos(2*pi*f1*t+1e-2*cumsum(usrDatRsm_1));


%==============================================%
% Quadrature oscillator

LO_i=cos(2*pi*f1*t);
LO_q=sin(2*pi*f1*t);

sFM_1_i=sFM_1.*LO_i;
sFM_1_q=sFM_1.*LO_q;

LPF=fir1(50,100e6/(Fd/2)); % If filter to remove the upconverted component
sFM_1_i=filter(LPF,1,sFM_1_i);
sFM_1_q=filter(LPF,1,sFM_1_q);

sFM_1_dem=atan2d(sFM_1_i,sFM_1_q);
% sFM_1_dem=unwrap(sFM_1_dem);
sFM_1_dem_diff=diff(sFM_1_dem);

%==============================================%
% Spectra 

[spectr_usrDatRsm_1, fr]=win_fft(usrDatRsm_1, 4e9,10^4,10^3);
[spectr_sFM_1_i, fr]=win_fft(sFM_1_i, 4e9,10^4,10^3);
[spectr_sFM_1_q, fr]=win_fft(sFM_1_q, 4e9,10^4,10^3);
[spectr_sFM_1_dem, fr]=win_fft(sFM_1_dem, 4e9,10^4,10^3);


%==============================================%
%% Plots

figure(1)
hold on
plot(t*1e6,usrDatRsm_1)
plot(t(1:end-1)*1e6,sFM_1_dem_diff)
% plot(t*1e6,sFM_1_dem)
ylim([-2, 2])
xlabel('t, ns')
ylabel('s(t), V')
legend("s_{1}(t)","s_{1}(t) Demodulated",'location','northeastOutside')
grid on, grid minor
set(gca,'fontsize',20)

figure(2)
hold on
plot(fr*1e-6,20*log10(spectr_sFM_1_i))
plot(fr*1e-6,20*log10(spectr_sFM_1_q))
plot(fr*1e-6,20*log10(spectr_usrDatRsm_1))
plot(fr*1e-6,20*log10(spectr_sFM_1_dem))
xlim([0, 500])
xlabel('f, MHz')
ylabel('s(f), dB')
legend("FM_{s(t)} I","FM_{s(t)} Q","s(t)","s(t) Demodulated",'location','northeast')
grid on, grid minor
set(gca,'fontsize',20)

