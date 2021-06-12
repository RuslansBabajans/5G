% Practical work 3
clear all; clc; close all

Fs=100e6; % baseband clock
N=50000; % number of symbols to transmit

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
gain_dist=2*ones(size(usrDatFlt));
% gain_dist(1:40000)=0;
usrRx=usrDatFlt.*gain_dist;
%-----------------------------------------

% Reference level calculations:
% E[(x^2-PW)*x^2]=E[x^4 -PWx^2]=E[x^4] -E[PWx^2]=E[x^4] -PW*E[x^2]=0     =>   PW=E[x^4]/E[x^2]
% E[x^2-PW]=E[x^2]-PW=0     =>   PW=E[x^2]
% E[|x|-PW]=E[|x|]-PW=0     =>   PW=E[|x|]


% Jumping window averager
mu=5e-6;
% REF=mean(usrDatFlt.^4)/mean(usrDatFlt.^2);%4.6754;
% REF=mean(usrDatFlt.^2); % USE IN YOUR LABS
REF=mean(abs(usrDatFlt)); 
y1gain=zeros(size(usrRx));
y1data=zeros(size(usrRx));
gainAcc=1;
agcOut=0;
for k=1:length(usrRx)

    agcOut_tmp=agcOut;

    % Apply compensation
    agcOut=usrRx(k)*gainAcc;

    % Move parameter to the local minimum
    % gainAcc=gainAcc-8*mu*(agcOut_tmp^2-REF)*gainAcc*usrRx(k)^2;
    % gainAcc=gainAcc-8*mu*(agcOut_tmp^2-REF); % USE IN YOUR LABS
    gainAcc=gainAcc-8*mu*(abs(agcOut_tmp)-REF);

    % Output debug values
    y1gain(k)=gainAcc;
    y1data(k)=agcOut;


end




figure(1); hold on
plot(gain_dist,'b.-')
plot(y1gain,'g.-')
plot(gain_dist.*y1gain,'m.-')
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

