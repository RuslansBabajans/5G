close all, format compact, format long, clc,
%% Lab 1 Modulations and receivers

Fs = 100e6;     % Baseband clock
Fd = 4000e6;    % Analoge sampling freq
f0 = 170e6;    % Carrier freq
N = 200;        % Number of symbols

usrDat = kron(randi(2, 1, N)*2-3, ones(1, 4)); % symbol rate equal to 25Mbaud

% Resampling
usrDatRsm = resample(usrDat,40,1);

% Modulation
t = (1:length(usrDatRsm))/Fd;
sFM = cos(2*pi*f0*t+1e-2*cumsum(usrDatRsm));

% Calculate AM spectra
[spectr, fr] =     win_fft(usrDatRsm, 4e9,10^4,10^3);
[spectrFM, fr] =    win_fft(sFM, 4e9,10^4,10^3);

% 4dB attenuation
sFMatt = sFM; % .* 10^(-4/20);

% Level adjusting
P_FMin = 2*rms(usrDatRsm)^2;
P_FMatt = rms(sFMatt)^2;
sFMadj = sFMatt.*sqrt(P_FMin/P_FMatt);

% Mixer
I = cos(2*pi*f0*t);
sFMmxI = sFMadj .* I;
Q = sin(2*pi*f0*t);
sFMmxQ = sFMadj .* Q;
[spectrFMadj, fr]=win_fft(sFMadj, 4e9,10^4,10^3);
[spectrFMmxI, fr]=win_fft(sFMmxI, 4e9,10^4,10^3);
[spectrFMmxQ, fr]=win_fft(sFMmxQ, 4e9,10^4,10^3);

% baseband filter
LPF = fir1(70, 100e6/(Fd/2));
% fvtool(LPF)
sFMbbI = filter(LPF, 1, sFMmxI);
sFMbbQ = filter(LPF, 1, sFMmxQ);

[spectrFMbbI, fr]=win_fft(sFMbbI, 4e9,10^4,10^3);
[spectrFMbbQ, fr]=win_fft(sFMbbQ, 4e9,10^4,10^3);

% Instant Amp and phase
sFMamp = sqrt(sFMbbI.^2 + sFMbbQ.^2);
sFMphase = atan2d(sFMbbI, sFMbbQ);

[spectrFMamp, fr]=win_fft(sFMamp, 4e9,10^4,10^3);
[spectrFMphase, fr]=win_fft(sFMphase, 4e9,10^4,10^3);

% Derivative 
% sFMphase = rad2deg(sFMphase);
sFMdiff = diff(sFMphase);

% LPF
sFMdem = filter(LPF, 1, sFMdiff);

[spectrFMdem, fr]=win_fft(sFMdem, 4e9,10^4,10^3);


figure(), hold on, grid on, grid minor
plot(sFMdem, 'b.-')
plot(usrDatRsm, 'r.-')
hold off
return
figure(), hold on, grid on, grid minor;
plot(fr*1e-6, 20*log10(spectrFMmxQ), 'b.-'),
plot(fr*1e-6, 20*log10(spectrFMmxI), 'r.-'),
plot(fr*1e-6, 20*log10(spectrFMphase), 'g.-'),
plot(fr*1e-6, 20*log10(spectrFMdem), 'c.-'),
hold off,
return



figure(), hold on, grid on, grid minor
plot(usrDatRsm, 'b.-')
plot(sFM1flt, 'r.-')
hold off

figure(), hold on, grid on, grid minor
plot(usrDatRsm2, 'b.-')
plot(sFM2flt, 'r.-')
hold off
return

figure(), hold on, grid on, grid minor
plot(fr*1e-6, 20*log10(spectrFMdem1), 'b.-')
plot(fr*1e-6, 20*log10(spectrFM1flt), 'r.-')
hold off

figure(), hold on, grid on, grid minor
plot(fr*1e-6, 20*log10(spectrFMdem2), 'b.-')
plot(fr*1e-6, 20*log10(spectrFM2flt), 'r.-')
hold off

return