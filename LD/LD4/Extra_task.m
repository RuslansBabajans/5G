% Practical work 3
clear all; clc; close all

Fs=100e6; % baseband clock
N=50000; % number of symbols to transmit

% K=5000; % aevraging window length

% Create alphabeth
ALPH=-3:2:3;
ALPH=repmat(-3:2:3,4,1)+1j*repmat(ALPH.',1,4);

% Pulse-shaping filter
firTx=firrcos(64,0.25,0.25,2,'rolloff','normal'); firTx=firTx/max(firTx); firRx=1;

%-------------- Transmitter --------------
% Generate user data
usrDat=kron(ALPH(randi(numel(ALPH),1,N)),[1 0 0 0]);

% Transmitter filtration - pulse shaping
usrDatFlt=filter(firTx, 1, usrDat); usrDatFlt=usrDatFlt(65:end);


% Crate phase offset signal and add into the used signal
phase_dist=pi/6*ones(size(usrDatFlt));
% gain_dist(1:40000)=0;
usrRx=usrDatFlt.*exp(j*phase_dist);
%-----------------------------------------

usrRx=usrRx(1:4:end);
usrDatFlt=usrDatFlt(1:4:end);
phase_dist=phase_dist(1:4:end);

% plot(usrRx,'b.')
% return

% Jumping window averager
mu=1e-3;
y1phase=zeros(size(usrRx));
y1data=zeros(size(usrRx));
phaseAcc=1;
phseOut=0;
phseEst=0;
for k=1:length(usrRx)

    phseOut_tmp=phseOut;
    phseEst_tmp=phseEst;

    % Apply compensation
    phseOut=usrRx(k)*exp(-j*phaseAcc);
    phseEst=usrDatFlt(k);

    % Move parameter to the local minimum
    phaseAcc=phaseAcc+real(2j*mu*(phseOut_tmp-phseEst_tmp)*phseOut_tmp);


    % Output debug values
    y1phase(k)=phaseAcc;
    y1data(k)=phseOut;


end


figure(1); hold on
plot(phase_dist,'b.-')
plot(y1phase,'g.-')
plot(phase_dist-y1phase,'m.-')
% plot(usrDatFlt-y1data,'r.-')

% return
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

