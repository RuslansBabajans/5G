%% Clear variables and close figures
format long
clear variables
close all
%==============================================%
%% Parameters
Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
f2=7e6; % to subtract
f1=10e6; % to add
%==============================================%
%% Generate user data
t=(1:1e6)/Fd; % time vector
I_t=1*cos(2*pi*f1*t)+0.71*cos(2*pi*f2*t);
Q_t=1*sin(2*pi*f1*t)-0.71*sin(2*pi*f2*t);
%==============================================%
%% Upsampling block
usrDatRsmI=resample(I_t,40,1);
usrDatRsmQ=resample(Q_t,40,1);
%==============================================%
%% Quadrature modulation
% QAM is selected as s_mod(t)=I*cos-Q*sin
sAM=I_t.*cos(2*pi*f0*t)-Q_t.*sin(2*pi*f0*t); % QAM signal
%==============================================%
%%  QAM demodulator
sQAMdemI= sAM.*cos(2*pi*f0*t);
sQAMdemQ=-sAM.*sin(2*pi*f0*t);

% Low-pass filter
LPF=fir1(100,f0/(Fd/2));
sAMfltI=2*filter(LPF,1,sQAMdemI);
sAMfltQ=2*filter(LPF,1,sQAMdemQ);
sAMfltI=sAMfltI(161:end);
sAMfltQ=sAMfltQ(161:end);
%==============================================%
%%  Downsampling block
usrDatRsm2I=resample(sAMfltI(11:end),1,40);
usrDatRsm2Q=resample(sAMfltQ(11:end),1,40);
%==============================================%
%% Calculate spectra
[spectr, fr]=win_fft(sAM, 4e9,2e4,10^3);
%==============================================%
%% Plots
figure(1)
plot(fr*1e-6, 20*log10(spectr),'Linewidth',2)
grid on, grid minor
xlabel("f, MHz")
ylabel("QAM_{out}(f), dB")
set(gca, 'Xlim', [100 180], 'XTick', 100:5:180, 'XTickLabel', 100:5:180)
set(gca, 'Ylim', [-10 1], 'YTick', -10:1:1, 'YTickLabel', -10:1:1)
set(gca, 'fontsize', 15)

figure(2)
subplot(2,1,1)
plot(t(975005:end)*1e6,usrDatRsm2I,'Linewidth',2)
xlabel("t, \mus")
ylabel("I_{RX}(t), V")
set(gca, 'Xlim', [245  245.1])
grid on, grid minor
set(gca, 'fontsize', 15)

subplot(2,1,2)
plot(t(975005:end)*1e6,usrDatRsm2Q,'Linewidth',2)
xlabel("t, \mus")
ylabel("Q_{RX}(t), V")
set(gca, 'Xlim', [245  245.1])
grid on, grid minor
set(gca, 'fontsize', 15)
