%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters

Fs=100e6; % 100 MHz baseband clock
Fd=4000e6; % analogue sampling freq
N=200; % number of symbols to transmit


% Generate user data
% usrDat=kron(randi(2,1,N)*2-3,ones(1,4));
usrDat=kron(randi(2,1,N)*2-3,ones(1,4));    % ones(1,4) is 4 samples per symbol,  Fs/4= 25 Mbaud

% Resampling
usrDatRsm=resample(usrDat,40,1);

% Time vector
t=(1:length(usrDatRsm))/Fd;
%==============================================%
%% Crystal receiver

f0=140e6; % Carrier frequency


% Modulation
sAM_TX=(2+usrDatRsm).*cos(2*pi*f0*t);

% Attenuation
sAM_RX=sAM_TX*0.630957; % 0.630957 corresponds to -4 dB voltage attenuation 
attenuation_check=20*log10(sAM_RX/sAM_TX); % attenuation_check= -4.0000047 dB

% AM demodulation using crystal receiver
sAMdem=abs(sAM);


