% Practical work 6 Equalization
clear all; clc; close all

Fs=100e6; % baseband clock
f0=1e9; % carrier frequency
N=10000; % number of symbols to transmit

% Create alphabeth
ALPH=-3:2:3; ALPH=repmat(-3:2:3,4,1)+1j*repmat(ALPH.',1,4);

% Pulse-shaping filter
firTx=firrcos(64,0.25,0.25,2,'rolloff','normal'); firTx=firTx/max(firTx);

%-------------- Transmitter --------------
% Generate user data
usrDat=kron(ALPH(randi(numel(ALPH),1,N)),[1 0 0 0]);

% Transmitter filtration - pulse shaping
usrDatFlt=filter(firTx, 1, usrDat); usrDatFlt=usrDatFlt(65:end);
%-----------------------------------------

% figure(1); hold on
% plot(usrDatFlt(1:4:end),'b.')

%% Two-ray propagation channel model
% tau=5ns;

DP=10^(-3/20);

% h=1*[zeros(1,8) 1 zeros(1,9)]+(1-DP)*lpntrp(17,0.5)*exp(-j*2*pi*f0*5.5e-9);
h=1/sqrt(1+(1-DP)^2)*[zeros(1,8) 1 zeros(1,9)]+(1-DP)/sqrt(1+(1-DP)^2)*lpntrp(17,0.5)*exp(-j*2*pi*f0*5.5e-9);
%fvtool(h)

usrChan=filter(h,1,usrDatFlt);
usrChan=usrChan(9:end);
usrChan=usrChan/sqrt(mean(abs(usrChan).^2))*sqrt(mean(abs(usrDatFlt).^2));

figure(1); hold on
plot(usrChan(1:4:end),'b.')
% return

% fvtool(usrChan,1,h,1)

%% Complex Gaussian process generation
% x=1*randn(1,N)+0+1j*randn(1,N)-0j;
% hist(abs(x),100)



%% Equalization

% Dispersion constant for CMA
R=mean(abs(ALPH).^4)/mean(abs(ALPH).^2); 

cntMain=0;
% eqCoef=zeros(1,15); eqCoef(9)=1;
eqCoef=zeros(1,33); eqCoef(17)=1;
eqFltr=zeros(size(eqCoef));
yOut=0;
yCst=0;

% Debug signals
yOutDat=zeros(size(usrChan));
yOutCnv=zeros(size(usrChan));

for k=1:length(usrChan)

    % Equalization
    yOut=eqFltr*eqCoef.';
    
    % Coefficient adjustment
    if cntMain==0
        yCst=yOut-ALPH(find(min(abs(yOut-ALPH(:).'))==abs(yOut-ALPH(:).'),1));
        eqCoef=eqCoef-0.001*yCst*conj(eqFltr); % LMS
        % eqCoef=eqCoef-0.00001*(abs(yOut)^2-R)*yOut*conj(eqFltr); % CMA
    end
      
    % Filter delay line
    eqFltr=[usrChan(k) eqFltr(1:end-1)];
        
    % Counter to indicate symbols positions
    cntMain=mod(cntMain+1,4);
    
    % Form output
    yOutDat(k)=yOut;
    yOutCnv(k)=yCst;
    
end


figure(2); hold on
subplot(2,2,1); plot(yOutDat(1+10000:4:end),'b.')
subplot(2,2,2); plot(yOutDat(2+10000:4:end),'b.')
subplot(2,2,3); plot(yOutDat(3+10000:4:end),'b.')
subplot(2,2,4); plot(yOutDat(4+10000:4:end),'b.')

figure(3); hold on
plot(20*log10(abs(yOutCnv(1:4:end))),'b.')

