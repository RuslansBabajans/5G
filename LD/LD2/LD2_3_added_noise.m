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
firTx=firrcos(46,0.5,0.25,2,'rolloff','normal');
usrDatFltI=filter(firTx, 1, usrDatI);
usrDatFltQ=filter(firTx, 1, usrDatQ);
usrDatFltI=usrDatFltI(46:end); % Skip transition time
usrDatFltQ=usrDatFltQ(46:end);
%==============================================%
%% Upsampling block
usrDatRsmI=resample(usrDatFltI,40,1); % Resample to 4 Ghz sample rate
usrDatRsmQ=resample(usrDatFltQ,40,1);
t=(1:length(usrDatRsmI))/Fd; % time vector
I_t=usrDatRsmI;
Q_t=usrDatRsmQ;
%==============================================%
%% Quadrature modulation
% QAM is selected as s_mod(t)=I*cos-Q*sin
sAM=I_t.*cos(2*pi*f0*t)-Q_t.*sin(2*pi*f0*t); % QAM signal
%==============================================%
%% AWGN channel model
SNR=20; % signal-to-noise ratio in dB
usrChan=sAM;
usrChan=usrChan+randn(size(usrChan))*sqrt(mean(abs(usrChan).^2)/2)*10^(-SNR/20);
%==============================================%
%%  QAM demodulator
sQAMdemI= usrChan.*cos(2*pi*f0*t);
sQAMdemQ=-usrChan.*sin(2*pi*f0*t);

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

fvtool(firTx)

figure(5)
plot(fr*1e-6,spectr_usrDatI,'Linewidth',2)
hold on
plot(fr*1e-6,spectr_usrDatFltI,'Linewidth',2)
grid on, grid minor
set(gca, 'Xlim', [0 100], 'XTick', 0:5:100, 'XTickLabel', 0:5:100)
xlabel("f, MHz")
ylabel("s(f), V")
legend("usrDat  before  raised-cosine  filter","usrDat  after  raised-cosine  filter")
set(gca, 'fontsize', 15)