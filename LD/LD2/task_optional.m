%% Clear variables and close figures
format long
clear variables
close all
%==============================================%
%% Parameters
Beta=0.5;   % Roll-off factor
Fc=0.25;    % Cutoff frequency
Ts=1/(2*Fc)*1e-5;    % Symbol's length 
f=0:1e3:1e5; % Frequency vector
%==============================================%
%% Filter using frequency response expression from lecture 4
f_1=f(f<(1-Beta)/(2*Ts));
f_2=f((f>=(1-Beta)/(2*Ts))&(f<(1+Beta)/(2*Ts)));
f_3=f(f>=(1+Beta)/(2*Ts)); 
 
H_1=ones(size(f_1));
H_2=1/2*(1+cos((pi*Ts/Beta).*(abs(f_2)-((1-Beta)/(2*Ts)))));
H_3=zeros(size(f_3));

H_f=[H_1,H_2,H_3];
%==============================================%
%% Filter using firrcos
firTx=firrcos(96,Fc,Beta,2,'rolloff', 'normal');
[H_filter_all w_filter_all]=freqz(firTx); % Full frequency response of the filter
f_filter_all=abs(w_filter_all/(pi)); 
H_filter_all=(abs(H_filter_all));
%==============================================%
%% Plots
figure(2)
plot(f_filter_all, H_filter_all, 'b.-')
hold on
plot(f/f(end),H_f,'r.--');
xlabel("Normalized frequency (\times \pi rad/sample)")
ylabel("Magnitude, V")
legend("'firrcos' produced filter","Filter from the frequency response expression")
set(gca, 'fontsize', 12)
set(gca, 'Xlim', [0 1], 'XTick', 0:0.1:1, 'XTickLabel', 0:0.1:1)
set(gca, 'Ylim', [0 1.1], 'YTick', 0:0.1:1.1, 'YTickLabel', 0:0.1:1.1)
grid on, grid minor