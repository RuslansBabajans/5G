% Practical work 3
clear all; clc; close all

Fs=100e6; % baseband clock
N=100000; % number of symbols to transmit

% K=5000; % aevraging window length

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
dc_dist=1*ones(size(usrDatFlt));
dc_dist(1:10000)=0;
usrRx=usrDatFlt+dc_dist;
%-----------------------------------------




% Jumping window averager
mu=1/50000;
% y1dc=zeros(size(usrRx));
y1data=zeros(size(usrRx));
avgAcc=0;
avgOut=0;
for k=1:length(usrRx)

    avgOut_tmp=avgOut;

    % Apply compensation
    avgOut=usrRx(k)-avgAcc;

    % Move parameter to the local minimum
    avgAcc=avgAcc+2*mu*avgOut_tmp;

    % Output debug values
    % y1dc(k)=avgOut/K;
    y1data(k)=avgOut;


end




figure(1); hold on
plot(dc_dist,'b.-')
plot(usrDatFlt-y1data,'r.-')

figure(2); hold on
plot(usrDatFlt(1:1000),'b.-')
plot(y1data(1:1000),'r.-')

figure(3); hold on
plot(usrDatFlt(end-1000:end),'b.-')
plot(y1data(end-1000:end),'r.-')



% figure(2); hold on
% plot(dc_dist,'b.-')
% plot(y1dc,'r.-')
% plot(y2dc,'g.-')
% plot(y3dc,'m.-')

