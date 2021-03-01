%% Clear variables and close figures
format long
clear variables
close all

%==============================================%
%% Parameters

Fs=100e6; % 100 MHz baseband clock

N=200; % number of symbols to transmit


% Generate user data
% usrDat=kron(randi(2,1,N)*2-3,ones(1,4));
usrDat=kron(randi(2,1,N)*2-3,ones(1,4));
