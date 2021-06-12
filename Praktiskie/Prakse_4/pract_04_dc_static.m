% Practical work 3
clear all; clc; close all

Fs=100e6; % baseband clock
N=50000; % number of symbols to transmit

K=5000; % aevraging window length

% Create alphabeth
ALPH=-3:2:3;



% Pulse-shaping filter
firTx=firrcos(64,0.25,0.25,2,'rolloff','normal'); firTx=firTx/max(firTx); firRx=1;

%-------------- Transmitter --------------
% Generate user data
usrDat=kron(ALPH(randi(numel(ALPH),1,N)),[1 0 0 0]);

% Transmitter filtration - pulse shaping
usrDatFlt=filter(firTx, 1, usrDat); usrDatFlt=usrDatFlt(64:end);


% Crate DC-offset signal and add it to the used signal
dc_dist=0.5*ones(size(usrDatFlt));
dc_dist(1:40000)=0;
usrRx=usrDatFlt+dc_dist;
%-----------------------------------------




% Jumping window averager
y1dc=zeros(size(usrRx));
y1data=zeros(size(usrRx));
avgCnt=0;
avgAcc=0;
avgOut=0;
for k=1:length(usrRx)


    % Output debug values
    y1data(k)=usrRx(k)-avgOut/K;
    y1dc(k)=avgOut/K;

    % Acummulator and its output
    if avgCnt==0
        avgOut=avgAcc;
        avgAcc=usrRx(k);
    else
        avgOut=avgOut;
        avgAcc=avgAcc+usrRx(k);
    end


    % Counter for window length
    avgCnt=mod(avgCnt+1,K);

end



% Sliding window averager
y2dc=zeros(size(usrRx));
y2data=zeros(size(usrRx));
avgMem=zeros(1,K);
avgAcc=0;
for k=1:length(usrRx)


    % Output debug values
    y2data(k)=usrRx(k)-avgAcc;
    y2dc(k)=avgAcc;

    % Acummulator and its output
    avgAcc=avgAcc+(usrRx(k)-avgMem(end))/K;

    % Memory to keep samples
    avgMem=[usrRx(k) avgMem(1:end-1)];

end


% IIR averager
y3dc=zeros(size(usrRx));
y3data=zeros(size(usrRx));
avgAcc=0;
for k=1:length(usrRx)


    % Output debug values
    y3data(k)=usrRx(k)-avgAcc;
    y3dc(k)=avgAcc;

    % Acummulator in IIR
    avgAcc=(1-1/K)*avgAcc+usrRx(k)/K;


end


% figure(1); hold on
% plot(usrDatFlt,'b.-')
% plot(y1data,'r.-')
% plot(y2data,'g.-')
% plot(y3data,'m.-')


figure(2); hold on
plot(dc_dist,'b.-')
plot(y1dc,'r.-')
plot(y2dc,'g.-')
plot(y3dc,'m.-')
