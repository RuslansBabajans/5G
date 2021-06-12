% Practical work 3
clear all; clc; close all

Fs=100e6; % baseband clock
% f0=140e6; % carrier frequency
N=30000; % number of symbols to transmit

% Create alphabeth
ALPH=-3:2:3;
ALPH=repmat(-3:2:3,4,1)+1j*repmat(ALPH.',1,4);

% Pulse-shaping filter
firTx=firrcos(64,0.25,0.25,2,'rolloff','normal'); firTx=firTx/max(firTx); firRx=1;
% firTx=firrcos(64,0.25,0.25,2,'rolloff', 'sqrt'); firTx=4*firTx/sum(firTx); firRx=firTx/sum(firTx); % firTx=firTx/sqrt(max(conv(firTx,firTx))); firRx=firTx;

SNR=140; % AWGN signal-to-noise ratio in dB

fc=0e3; % carrier offset in Hz
phc=0*pi/180; % carrier phase shift in radians
phNs=0*pi/180; % phase noise standard deviation in radians (not standard dBc - just for visualization) 

txImpPhs=10*pi/180; % Tx IQ-imbalance phase in radians -pi/12...pi/12
txImpDCi=0; % Tx IQ-imbalance I-component DC-offset -1...1
txImpDCq=0; % Tx IQ-imbalance Q-component DC-offset -1...1
txImpGni=1; % Tx IQ-imbalance I-component gain 0.5...1.5
txImpGnq=1; % Tx IQ-imbalance Q-component gain 0.5...1.5

BETA=0.0; % AM/AM nonlinearity parameter 0...1
GAMMA=-0.0; % AM/PM nonlinearity parameter -2...2

rxImpPhs=10*pi/180; % Rx IQ-imbalance phase in radians -pi/12...pi/12
rxImpDCi=0; % Rx IQ-imbalance I-component DC-offset -1...1
rxImpDCq=0; % Rx IQ-imbalance Q-component DC-offset -1...1
rxImpGni=1; % Rx IQ-imbalance I-component gain 0.5...1.5
rxImpGnq=1; % Rx IQ-imbalance Q-component gain 0.5..1.5


%-------------- Transmitter --------------
% Generate user data
usrDat=kron(ALPH(randi(numel(ALPH),1,N)),[1 0 0 0]);

% Transmitter filtration - pulse shaping
usrDatFlt=filter(firTx, 1, usrDat); usrDatFlt=usrDatFlt(33:end);
% fvtool(firTx)
% return
% figure(), hold on
% plot(abs(usrDat(1:250)),'b.-')
% plot(abs(usrDatFlt(1:250)),'r.-')
% plot(abs(usrDatFlt(64:314)),'g.-')
% return
%-----------------------------------------
% pre-distortion
txImpPhsComp = -txImpPhs;
usrTxImp1=txImpGni*real(usrDatFlt)+1j*txImpGnq*imag(usrDatFlt)+txImpDCi+1j*txImpDCq;
usrTxImp1=cos(txImpPhsComp)*real(usrTxImp1)+sin(txImpPhsComp)*imag(usrTxImp1)+1j*(sin(txImpPhsComp)*real(usrTxImp1)+cos(txImpPhsComp)*imag(usrTxImp1));

% Transmitter IQ-imbalance
usrTxImp=txImpGni*real(usrDatFlt)+1j*txImpGnq*imag(usrDatFlt)+txImpDCi+1j*txImpDCq;
usrTxImp=cos(txImpPhs)*real(usrTxImp)+sin(txImpPhs)*imag(usrTxImp)+1j*(sin(txImpPhs)*real(usrTxImp)+cos(txImpPhs)*imag(usrTxImp));
usrTxImp = (usrTxImp + usrTxImp1)./2;

% Nonlinear distortions
RAD=unique(abs(ALPH)); NORM=3.5*sqrt(2);
usrNonlin=usrTxImp./(1+BETA*abs(usrTxImp/NORM).^2).*exp(1j*(GAMMA*abs(usrTxImp/NORM).^2)./(1+abs(usrTxImp/NORM).^2));
usrNonlin=usrNonlin.*exp(-1j*GAMMA*abs(RAD(1)/NORM).^2/(1+abs(RAD(1)/NORM).^2)); % inner radius opposite rotation

% Carrier desynchronization in the demodulation process
t=(1:length(usrNonlin))/Fs;
usrChan=usrNonlin.*exp(j*2*pi*fc*t+j*phc+j*phNs*randn(1,length(t)));
usrChan=usrChan+(randn(size(usrChan))+1j*randn(size(usrChan)))*sqrt(mean(abs(usrChan).^2)/2)*10^(-SNR/20);

% Receiver IQ-imbalance
usrRxImp=cos(rxImpPhs)*real(usrChan)-sin(rxImpPhs)*imag(usrChan)+1j*(-sin(rxImpPhs)*real(usrChan)+cos(rxImpPhs)*imag(usrChan));
usrRxImp=rxImpGni*real(usrRxImp)+1j*rxImpGnq*imag(usrRxImp)+rxImpDCi+1j*rxImpDCq;

% post compensation
rxImpPhsComp = -rxImpPhs;
% usrRxImp2=cos(rxImpPhsComp)*real(usrChan)-sin(rxImpPhsComp)*imag(usrChan)+1j*(-sin(rxImpPhsComp)*real(usrChan)+cos(rxImpPhsComp)*imag(usrChan));
usrRxImp2=cos(rxImpPhsComp)*real(usrRxImp)-sin(rxImpPhsComp)*imag(usrRxImp)+1j*(-sin(rxImpPhsComp)*real(usrRxImp)+cos(rxImpPhsComp)*imag(usrRxImp));
% usrRxImp2=rxImpGni*real(usrRxImp2)+1j*rxImpGnq*imag(usrRxImp2)+rxImpDCi+1j*rxImpDCq;
usrRxImp =  usrRxImp2; %(usrRxImp + usrRxImp2)./2;

%---------------- Receiver ---------------
% Clipping
clpLev=4;
usrClip=usrRxImp;
% usrClip(abs(real(usrRxImp))>clpLev)=clpLev*sign(real(usrClip(abs(real(usrRxImp))>clpLev)))+1j*imag(usrClip(abs(real(usrRxImp))>clpLev));
% usrClip(abs(imag(usrRxImp))>clpLev)=1j*clpLev*sign(imag(usrClip(abs(imag(usrRxImp))>clpLev)))+real(usrClip(abs(imag(usrRxImp))>clpLev));

% Matched filtration
usrRxFltr=filter(firRx, 1, usrClip); usrRxFltr=usrRxFltr(64:end);

% Automatic gain control
% usrRxAgc=usrRxFltr;
usrRxAgc=usrRxFltr*sqrt(mean(abs(usrDatFlt).^2)/mean(abs(usrRxFltr).^2));

% Carrier recovery
% usrDeRot=usrRxAgc;
usrDeRot=usrRxAgc.*exp(-j*2*pi*fc*t(64:end)-j*phc);
%-----------------------------------------

% eyediagram(usrDatRsm2Q(3:end),16)
% return
a = 2;
figure(); axis square; hold on; grid on
plot(usrDeRot(a:4:end),'b.')
plot(ALPH,'ro')
xlim([-4 4]); ylim([-4 4]),
xlabel('In-phase'), ylabel('Quadrature')
% return

% Calculate AM spectra
[spectr, fr]=win_fft(usrDatFlt(64:10063), 2,1000,10);
[spectr1, fr]=win_fft(usrDeRot(1:10000), 2,1000,10);
% [spectrAM, fr]=win_fft(sAM, 4e9,10^4,10^3);
%[spectrQAMdem, fr]=win_fft(sQAMdemQ, 4e9,10^4,10^3);
% [spectrQAMflt, fr]=win_fft(sAMfltQ, 4e9,10^4,10^3);


figure(); hold on; grid on
plot(fr, 20*log10(spectr),'b.-')
plot(fr, 20*log10(spectr1),'r.-')
% plot(fr, 20*log10(spectrQAMdem),'g.-')
% plot(fr, 20*log10(spectrQAMflt),'m.-')
% return

figure(), hold on
subplot(2,1,1), hold on
plot(real(usrDatFlt(64:1063)),'b.-')
plot(real(usrDeRot(1:1000)),'r.-')
title("In-phase")
hold off

subplot(2,1,2), hold on
plot(imag(usrDatFlt(64:1063)),'b.-')
plot(imag(usrDeRot(1:1000)),'r.-')
title("Quadrature")
hold off
