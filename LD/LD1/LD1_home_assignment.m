%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters
Fs = 1e3; % Sampling frequency                    
tau = 1/Fs; % Sampling period

m=0.5; % modulation depth
Sm=1; % Carrier amplitude
delta_Sm=Sm*m;

f_c=Fs/20; % carrier frequency in Hz
f_m=f_c/10; % modulating signal frequency in Hz

t = 0:tau:2*(1/f_m);% Time vector

%==============================================%
% Signals
s_mod=delta_Sm*cos(2*pi*f_m.*t);
s_c=Sm.*sin(2*pi*f_c.*t);
s_AM_t=Sm.*(1+m*cos(2*pi*f_m.*t)).*cos(2*pi*f_c.*t); % AM signal

%==============================================%
%% Time plots
figure(1)
subplot(3,1,1)
plot(t,s_mod,'Linewidth',2)
xlabel('t, s')
ylabel('s_{mod}(t), V')
grid on, grid minor
set(gca,'fontsize',12)

subplot(3,1,2)
plot(t,s_c,'Linewidth',2)
xlabel('t, s')
ylabel('s_{c}(t), V')
grid on, grid minor
set(gca,'fontsize',12)

subplot(3,1,3)
plot(t,s_AM_t,'Linewidth',2)
xlabel('t, s')
ylabel('s_{AM}(t), V')
grid on, grid minor
set(gca,'fontsize',12)
%==============================================%
%% Spectrum plots
f=[f_c-f_m      f_c     f_c+f_m];
s_AM_f=[0.5*m*Sm     Sm     0.5*m*Sm];

figure(2)
stem(f,s_AM_f,'Linewidth',2)
hold on
stem(f_m,delta_Sm,'Linewidth',2)

xlabel('f, Hz')
ylabel('s(f), V')
grid on, grid minor
set(gca,'fontsize',12)
ylim([0, 1.1])
set(gca, 'XLim', [0, 60], 'XTick', 0:5:60,...
    'XTickLabel', 0:5:60);
legend("s_{AM}(f)","s_{mod}(f)",'location','northeast')



