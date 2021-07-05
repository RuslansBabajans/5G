% LD 6 task 1
%% Clear variables and close figures
format long
clear variables
close all
%=========================================================================%
%% Parameters
Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=1e9; % carrier frequency
N=1000; % number of symbols to transmit
ALPH=-3:2:3;
%=========================================================================%
%% Tx side
    %% Generate user data
    % 5 samples per symbol, Fs/4= 25 Mbaud
    usrDatI=kron(ALPH(randi(4,1,N)),[1 0 0 0]);  
    usrDatQ=kron(ALPH(randi(4,1,N)),[1 0 0 0]);
    %=====================================================================%
    %% Raised-cosine filters 
    firTx=firrcos(64,0.25,0.25,2,'rolloff','normal');
    firTx=firTx/max(firTx);
    firRx=1;
    %=====================================================================%
    %% Signal forming using a single raised-cosine filter
    usrDatFltI=filter(firTx, 1, usrDatI);
    usrDatFltQ=filter(firTx, 1, usrDatQ);

    usrDatFltI=usrDatFltI(66:end); % Skip filter transition time
    usrDatFltQ=usrDatFltQ(66:end);
    %=====================================================================%
    %% Upsampling block
    usrDatRsmI=resample(usrDatFltI,40,1); % Resample to 4 Ghz sample rate
    usrDatRsmQ=resample(usrDatFltQ,40,1);
    t=(1:length(usrDatRsmI))/Fd; % time vector
    I_t=usrDatRsmI;
    Q_t=usrDatRsmQ;
    %=====================================================================%
    %% Quadrature modulation
    % QAM is selected as s_mod(t)=I*cos-Q*sin
    sAM=I_t.*cos(2*pi*f0*t)-Q_t.*sin(2*pi*f0*t); % QAM signal
%=========================================================================%
%% Channel models
    %% Direct AWGN channel with no time delay
    SNR=120; % Signal-to-noise ratio in dB
    usrChan=sAM;
    usrChan=usrChan+randn(size(usrChan))*sqrt(mean(abs(usrChan).^2)/2)...
            *10^(-SNR/20);
    %=====================================================================%
    %% Two-ray propagation channel model
    DP=10^(-5/20);
    
%     h=1*[zeros(1,8) 1 zeros(1,9)]+(1-DP)*lpntrp(17,0.5)...
%       *exp(-j*2*pi*f0*5.5e-9);
    
    h=1/sqrt(1+(1-DP)^2)*[zeros(1,8) 1 zeros(1,9)]+...
      (1-DP)/sqrt(1+(1-DP)^2)*lpntrp(17,0.5)*exp(-j*2*pi*f0*5.5e-9);
  
  
    fvtool(h)
    %=====================================================================%
    
    
%=========================================================================%
%% Rx side
    %% QAM demodulator
    % Mixer
    LO=exp(-j*(2*pi*f0*t));   
    
    sQAMdemI=real(usrChan.*LO);
    sQAMdemQ=imag(-usrChan.*LO);
    % Low-pass filter
    LPF=fir1(100,f0/(Fd/2));
    
    sAMfltI=2*filter(LPF,1,sQAMdemI);
    sAMfltQ=2*filter(LPF,1,sQAMdemQ);
    
    sAMfltI=sAMfltI(161:end); % Skip filter transition time 
    sAMfltQ=sAMfltQ(161:end);      
    %=====================================================================%
    %% Downsampling block
    usrDatRsm2I=resample(sAMfltI(11:end),1,40);
    usrDatRsm2Q=resample(sAMfltQ(11:end),1,40);
    %=====================================================================%
    %% Symbol detection
    symbolI=usrDatRsm2I(5:4:end);
    symbolQ=usrDatRsm2Q(5:4:end);
%=========================================================================%
%% Spectra


%=========================================================================%
%% Plots  
figure(1)
plot(symbolI,symbolQ,'b.')
grid on, grid minor
xlabel("I")
ylabel("Q")


%=========================================================================%

