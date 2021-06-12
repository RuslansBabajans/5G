%% Clear variables and close figures
format long
clear variables
close all
%==============================================%
%% Parameters
Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
N=1000; % number of symbols to transmit
ALPH=-3:2:3;
%==============================================%
%% Generate user data
usrDatI=kron(ALPH(randi(4,1,N)),[1 0 0 0 0]);  % 5 samples per symbol, Fs/5= 20 Mbaud
usrDatQ=kron(ALPH(randi(4,1,N)),[1 0 0 0 0]);
%==============================================%
%% Pulse-shaping filter
firTx=firrcos(66,0.2,0.14,2,'rolloff', 'sqrt');

f_m=[0, 0.2, 0.25, 0.4, 0.7, 1];
H_m=[0, 0, -30, -40, -50, -50];
w_m=f_m*pi;

H_m_filter = freqz(firTx,1,w_m);
H_m_filter=round(20*log10(abs(H_m_filter)))

[H_filter_all w_filter_all]=freqz(firTx);
f_filter_all=w_filter_all/(pi);
H_filter_all=20*log10(abs(H_filter_all));

usrDatFltI=filter(firTx, 1, usrDatI);
usrDatFltQ=filter(firTx, 1, usrDatQ);

usrDatFltI=usrDatFltI(66:end); % Skip transition time
usrDatFltQ=usrDatFltQ(66:end);
%==============================================%
%% Upsampling block
usrDatRsmI=resample(usrDatFltI,40,1); % Resample to 4 Ghz sample rate
usrDatRsmQ=resample(usrDatFltQ,40,1);
t=(1:length(usrDatRsmI))/Fd; % time vector
I_t=usrDatRsmI;
Q_t=usrDatRsmQ;
%==============================================%
%% Quadrature modulation
phi_0=0*pi/180; 
% QAM is selected as s_mod(t)=I*cos-Q*sin
sAM=I_t.*cos(2*pi*f0*t+phi_0)-Q_t.*sin(2*pi*f0*t); % QAM signal
%==============================================%
%% AWGN channel model
SNR=120; % signal-to-noise ratio in dB
usrChan=sAM;
usrChan=usrChan+randn(size(usrChan))*sqrt(mean(abs(usrChan).^2)/2)*10^(-SNR/20);
%==============================================%
%% QAM demodulator
phi_1=15*pi/180; % pi/3
delta_f=0e3; % 100e3
f1=f0+delta_f;

LO_I=cos(2*pi*f1*t+phi_1);
LO_Q=sin(2*pi*f1*t);
sQAMdemI= usrChan.*LO_I;
sQAMdemQ=-usrChan.*LO_Q;

% -15         -0.324047031779238
% -10         -0.182014018891678
% -5          -0.078449071642888
%  0           0.026774724864121
%  5           0.158643876220417
% 10           0.222003277558050
% 15           0.247195672013946


% Jumping window averager
% mu=1e-3;
% y1phase=zeros(size(usrRx));
% y1data=zeros(size(usrRx));
% phaseAcc=1;
% phseOut=0;
% phseEst=0;
% 
% for k=1:length(usrRx)
% 
%     phseOut_tmp=phseOut;
%     phseEst_tmp=phseEst;
% 
%     % Apply compensation
%     phseOut=usrRx(k)*exp(-j*phaseAcc);
%     phseEst=usrDatFlt(k);
% 
%     % Move parameter to the local minimum
%     phaseAcc=phaseAcc+real(2j*mu*(phseOut_tmp-phseEst_tmp)*phseOut_tmp);
% 
% 
%     % Output debug values
%     y1phase(k)=phaseAcc;
%     y1data(k)=phseOut;
% 
% 
% end


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

usrDatRsm2I=filter(firTx, 1, usrDatRsm2I); usrDatRsm2I=usrDatRsm2I(60:end);
usrDatRsm2Q=filter(firTx, 1, usrDatRsm2Q); usrDatRsm2Q=usrDatRsm2Q(60:end);
%==============================================%
%% Calculate spectra
[spectr_usrDatI, fr]=win_fft(resample(usrDatI,40,1), 4e9,1e4,1e3);
[spectr_usrDatFltI, fr]=win_fft(resample(usrDatFltI,40,1), 4e9,1e4,1e3);
[spectr_usrChan, fr]=win_fft(usrChan, 4e9,1e4,1e3);
%==============================================%
%% Plots
figure(1)
plot(fr*1e-6, 20*log10(spectr_usrChan),'Linewidth',2)
grid on, grid minor
xlabel("f, MHz")
ylabel("QAM_{out}(f), dB")
set(gca, 'Xlim', [0 400], 'XTick', 0:20:400, 'XTickLabel', 0:20:400)
set(gca, 'fontsize', 15)
 
figure(2)
plot(t*1e6,sAM,'o-','Linewidth',2)
grid on, grid minor
xlabel("t, us")
ylabel("QAM_{out}(t), V")
set(gca, 'Xlim', [0 1], 'XTick', 0:0.1:1, 'XTickLabel', 0:0.1:1)
set(gca, 'fontsize', 15)

figure(3)
plot(usrDatRsm2I(5:5:end),usrDatRsm2Q(5:5:end),'b.')
grid on, grid minor
xlabel("I")
ylabel("Q")