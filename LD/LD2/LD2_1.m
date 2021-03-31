%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters

Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
% N=1000; % number of symbols to transmit
% 
% ALPH=-3:2:3;
f=5e6;
%==============================================%
%% Generate user data
t=(1:1e6)/Fd; % time vector

I_t=cos(2*pi*f*t);
Q_t=cos(2*pi*f*t);
% usrDatI=kron(ALPH(randi(4,1,N)),[1 0 0 0 0]);  % 5 samples per symbol, Fs/5= 20 Mbaud
% usrDatQ=kron(ALPH(randi(4,1,N)),[1 0 0 0 0]);

%==============================================%
%% Pulse-shaping filter

% firTx=firrcos(36,0.25,0.25,2,'rolloff','normal');
% 
% usrDatFltI=filter(firTx, 1, usrDatI);
% usrDatFltQ=filter(firTx, 1, usrDatQ);
% 
% usrDatFltI=usrDatFltI(36:end);
% usrDatFltQ=usrDatFltQ(36:end);

%==============================================%
%% Upsampling block

% usrDatRsmI=resample(usrDatFltI,40,1); % Resample to 4 Ghz sample rate
% usrDatRsmQ=resample(usrDatFltQ,40,1);

usrDatRsmI=resample(I_t,40,1);
usrDatRsmQ=resample(Q_t,40,1);
%==============================================%
%% Quadrature modulation
% QAM is selected as s_mod(t)=I*cos-Q*sin

% t=(1:length(usrDatRsmI))/Fd; % time vector
sAM=I_t.*cos(2*pi*f0*t)-Q_t.*sin(2*pi*f0*t); % QAM signal

%==============================================%
%%  QAM demodulator
% sQAMdemI= sAM.*cos(2*pi*f0*t);
% sQAMdemQ=-sAM.*sin(2*pi*f0*t);

% Low-pass filter

% LPF=fir1(100,f0/(Fd/2));
% sAMfltI=2*filter(LPF,1,sQAMdemI);
% sAMfltQ=2*filter(LPF,1,sQAMdemQ);
% sAMfltI=sAMfltI(161:end);
% sAMfltQ=sAMfltQ(161:end);
%==============================================%
%%  Downsampling block

% usrDatRsm2I=resample(sAMfltI(11:end),1,40);
% usrDatRsm2Q=resample(sAMfltQ(11:end),1,40);

%==============================================%
%% Calculate spectra

[spectr, fr]=win_fft(sAM, 4e9,10^4,10^3);

figure(1)
plot(fr*1e-6, 20*log10(spectr),'Linewidth',2)
grid on, grid minor
xlabel("f, MHz")
ylabel("QAM_{out}(f), dB")
set(gca, 'Xlim', [100 180], 'XTick', 100:5:180, 'XTickLabel', 100:5:180)
set(gca, 'fontsize', 15)


