%% Script to edit figures for LD3
% clear all; clc; close all
%====================================================%
figure(1)
grid minor
set(gca, 'fontsize', 12)
ylim([-4, 5])
legend("Compensated distortion","Origianal signal");
%====================================================%
figure(2)
grid minor
set(gca, 'fontsize', 12)
legend("Origianal signal","Compensated distortion");
xlabel("f, Hz")
ylabel("s(f), dB")
%====================================================%
figure(3)
subplot(2,1,1)
ylim([-6, 6])
grid minor
set(gca, 'fontsize', 12)
% legend("Origianal signal","Distorted signal");
xlabel("t, ms")
ylabel("Re[s(t)], V")


subplot(2,1,2)
ylim([-6, 6])
grid minor
set(gca, 'fontsize', 12)
% legend("Origianal signal","Distorted signal");
xlabel("t, ms")
ylabel("Im[s(t)], V")

%====================================================%
