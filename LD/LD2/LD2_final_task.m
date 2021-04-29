%% Clear variables and close figures
format long
% clear variables
close all
%==============================================%
%% Parameters
Fs=100e6; % baseband clock
Fd=4000e6; % analog sampling freq
f0=140e6; % carrier frequency
N=1000; % number of symbols to transmit
%==============================================%
%% Constellation points forming initials "RB"
ALPH_I=[ -1 -0.4 0.4 0.7,...
         -1 -0.5 0.4 1,...
         -1 -0.7 0.4 0.7,...
         -1 -0.7 0.4 0.7,...
         -1 -0.4 0.4 1,...
         -0.4  1,...
         -1 -0.7 0.4 0.7]*5;
     
ALPH_Q=[ -0.9 -0.9 -0.9 -0.9,...
         -0.6 -0.6 -0.6 -0.6,...
         -0.3 -0.3 -0.3 -0.3,...
            0    0    0    0,...
          0.3  0.3 0.3 0.3,...
          0.5  0.5,...
          0.6  0.6  0.6  0.6]*5;
%==============================================%
%% Form message signal
% Message signal= (500 random lowercase letters)+(fancy saying in latin)+(500 random lowercase letters)
% Form 1000 random letters split in 500+500
random_lowercase_letters_1=randi([97,122],1,250); % Lowercase ASCII letters
random_upperrcase_letters_1=randi([65,90],1,250); % Uppercase ASCII letters

random_misc_ascii_symbols_1=randi([32,64],1,100); % ASCII symbols
random_more_misc_ascii_symbols_1=randi([91,96],1,100); % ASCII symbols
random_even_more_misc_ascii_symbols_1=randi([123,126],1,50); % ASCII symbols

random_lowercase_letters_2=randi([97,122],1,250); % Lowercase ASCII letters
random_upperrcase_letters_2=randi([65,90],1,250); % Uppercase ASCII letters

random_misc_ascii_symbols_2=randi([32,64],1,100); % ASCII symbols
random_more_misc_ascii_symbols_2=randi([91,96],1,100); % ASCII symbols
random_even_more_misc_ascii_symbols_2=randi([123,126],1,50); % ASCII symbols

random_ascii_symbols_1_original=[random_lowercase_letters_1,...
    random_upperrcase_letters_1,...
    random_misc_ascii_symbols_1,...
    random_more_misc_ascii_symbols_1,...
    random_even_more_misc_ascii_symbols_1];

random_ascii_symbols_2_original=[random_lowercase_letters_2,...
    random_upperrcase_letters_2,...
    random_misc_ascii_symbols_2,...
    random_more_misc_ascii_symbols_2,...
    random_even_more_misc_ascii_symbols_2];

% Then edit generated letters and symbols: convert uppercase letters to
% lowercase, cut all other symbols (including spaces commas, numbers)

random_ascii_symbols_1_edited=random_ascii_symbols_1_original;
random_ascii_symbols_2_edited=random_ascii_symbols_2_original;

for n=length(random_ascii_symbols_1_original):-1:1
    
    if random_ascii_symbols_1_edited(n) >= 65 & random_ascii_symbols_1_edited(n) <= 90
       random_ascii_symbols_1_edited(n)=random_ascii_symbols_1_edited(n)+32; % Convert uppercase to lowercase
    end
    
    if (random_ascii_symbols_1_edited(n) < 65 | random_ascii_symbols_1_edited(n) > 122)...
      |(random_ascii_symbols_1_edited(n) > 90 & random_ascii_symbols_1_edited(n) < 97)
        random_ascii_symbols_1_edited(n)=[]; % Cut ASCII symbol
    end
    
    %--------------------------------------------------------------------
    
    if random_ascii_symbols_2_edited(n) >= 65 & random_ascii_symbols_2_edited(n) <= 90
       random_ascii_symbols_2_edited(n)=random_ascii_symbols_2_edited(n)+32; % Convert uppercase to lowercase
    end
    
    if (random_ascii_symbols_2_edited(n) < 65 | random_ascii_symbols_2_edited(n) > 122)...
      |(random_ascii_symbols_2_edited(n) > 90 & random_ascii_symbols_2_edited(n) < 97)
        random_ascii_symbols_2_edited(n)=[]; % Cut ASCII symbol  
    end
end

% char(random_ascii_symbols_1_original)
% char(random_ascii_symbols_1_edited)
% length(random_ascii_symbols_1_edited)

% Form the fancy message and process same as the random letter sequencies
fancy_phrase_in_latin=char('FECI QUOD POTUI, FACIANT MELIORA POTENTESA');
fancy_phrase_in_latin_to_decimal=double(fancy_phrase_in_latin);
fancy_phrase_in_latin_decimal_edited=fancy_phrase_in_latin_to_decimal;

for n=length(fancy_phrase_in_latin_to_decimal):-1:1
    if fancy_phrase_in_latin_decimal_edited(n) >= 65 & fancy_phrase_in_latin_decimal_edited(n) <= 90
        fancy_phrase_in_latin_decimal_edited(n)=fancy_phrase_in_latin_decimal_edited(n)+32; % Convert uppercase to lowercase
    end
    
    if (fancy_phrase_in_latin_decimal_edited(n) < 65 | fancy_phrase_in_latin_decimal_edited(n) > 122)...
      |(fancy_phrase_in_latin_decimal_edited(n) > 90 & fancy_phrase_in_latin_decimal_edited(n) < 97)
        fancy_phrase_in_latin_decimal_edited(n)=[]; % Cut ASCII symbol 
    end
end
 
message_signal_ascii=[random_ascii_symbols_1_edited,...
                      fancy_phrase_in_latin_decimal_edited,...
                      random_ascii_symbols_2_edited]; % Form full message

message_signal_for_symbols=message_signal_ascii-96; % Convert ASCII decimal values to be represented with QAM symbols  

Data_string_I=ALPH_I(message_signal_for_symbols); % Convert to QAM In-phase values
Data_string_Q=ALPH_Q(message_signal_for_symbols); % Convert to QAM Quadrature values
%==============================================%
%% Convert user data to 20 MBd symbol rate
usrDatI=kron(Data_string_I,[1 0 0 0 0]); % 5 samples per symbol, Fs/5= 20 Mbaud
usrDatQ=kron(Data_string_Q,[1 0 0 0 0]);
%==============================================%
%% Pulse-shaping filter
% Create  Root-raised-cosine filter that complies with the spectral mask
firTx=firrcos(66,0.2,0.14,2,'rolloff', 'sqrt');

f_m=[0, 0.2, 0.25, 0.4, 0.7, 1]; % Spectral mask frequency
H_m=[0, 0, -30, -40, -50, -50]; % Spectral mask impulse response value
w_m=f_m*pi; % convert to Nyquist rad/sample; 

H_m_filter = freqz(firTx,1,w_m); % Filter response samples at w_m points
H_m_filter=round(20*log10(abs(H_m_filter))) % Convert to dB

[H_filter_all w_filter_all]=freqz(firTx); % Full frequency response of the filter 512 points by default
f_filter_all=w_filter_all/(pi); 
H_filter_all=20*log10(abs(H_filter_all));

% Tx Root-raised-cosine filter
usrDatFltI=filter(firTx, 1, usrDatI);
usrDatFltQ=filter(firTx, 1, usrDatQ);

usrDatFltI=usrDatFltI(66:end); % Skip transition time
usrDatFltQ=usrDatFltQ(66:end);
%==============================================%
%% Upsampling block
usrDatRsmI=resample(usrDatFltI,40,1); % Resample to 4 Ghz sample rate
usrDatRsmQ=resample(usrDatFltQ,40,1);
t=(1:length(usrDatRsmI))/Fd; % time vector
%==============================================%
%% Quadrature modulation
% QAM is selected as s_mod(t)=I_t*cos-Q_t*sin
I_t=usrDatRsmI;
Q_t=usrDatRsmQ;
sAM=(I_t+j.*Q_t)*0.5.*exp(j*2*pi*f0*t)+(I_t-j*Q_t)*0.5.*exp(-j*2*pi*f0*t);
%==============================================%
%% AWGN channel model
SNR=120; % signal-to-noise ratio in dB
usrChan=sAM;
usrChan=usrChan+randn(size(usrChan))*sqrt(mean(abs(usrChan).^2)/2)*10^(-SNR/20);
%==============================================%
%% QAM demodulator
% Local oscillator model
phi_0=0; % pi/3
delta_f=0; % 100e3
f1=f0+delta_f;
LO=exp(j*(2*pi*f1*t+phi_0));

% Demodulation
sQAMdemI=real(usrChan.*LO);
sQAMdemQ=imag(-usrChan.*LO);

% Low-pass filter
LPF=fir1(100,f0/(Fd/2));
sAMfltI=2*filter(LPF,1,sQAMdemI);
sAMfltQ=2*filter(LPF,1,sQAMdemQ);
sAMfltI=sAMfltI(161:end);
sAMfltQ=sAMfltQ(161:end);
%==============================================%
%%  Downsampling block
usrDatRsm2I=resample(sAMfltI(11:end),1,40);
usrDatRsm2Q=resample(sAMfltQ(11:end),1,40);

% Rx Root-raised-cosine filter
usrDatRsm2I=filter(firTx, 1, usrDatRsm2I);
usrDatRsm2I=usrDatRsm2I(60:end);
usrDatRsm2Q=filter(firTx, 1, usrDatRsm2Q);
usrDatRsm2Q=usrDatRsm2Q(60:end);
%==============================================%
%%  Detect message signal
% Round the symbols
usrDat_Rx_I_round=round(usrDatRsm2I, 1);
usrDat_Rx_Q_round=round(usrDatRsm2Q, 1);

% The actual symbol is every 5th sample
usrDat_Rx_I_symbols=usrDat_Rx_I_round(5:5:end);
usrDat_Rx_Q_symbols=usrDat_Rx_Q_round(5:5:end);

% Convert symbols to ASCII decimal
Message_signal_Rx_ascii_decimals=zeros(1,length(usrDat_Rx_I_symbols));
for n=1:1:length(usrDat_Rx_I_symbols)
    for k=1:1:length(ALPH_I)
        if usrDat_Rx_I_symbols(n) == ALPH_I(k)/5 &...
           usrDat_Rx_Q_symbols(n) == ALPH_Q(k)/5
       
            Message_signal_Rx_ascii_decimals(n)=k;
        end
    end
end

% Convert ASCII decimals to chars
Message_signal_Rx_ascii_chars=char(Message_signal_Rx_ascii_decimals+96);

% Find and display the fancy latin saying
% "feciquodpotuifaciantmeliorapotentesa" in the received message
latin_message_start_point=strfind(Message_signal_Rx_ascii_chars,...
                          'feciquodpotuifaciantmeliorapotentesa');
Recovered_message=Message_signal_Rx_ascii_chars(latin_message_start_point:latin_message_start_point+35)
%==============================================%
%% Calculate spectra
[spectr_usrDatI, fr]=win_fft(resample(usrDatI,40,1), 4e9,1e4,1e3);
[spectr_usrDatFltI, fr]=win_fft(resample(usrDatFltI,40,1), 4e9,1e4,1e3);
[spectr_usrChan, fr]=win_fft(usrChan, 4e9,1e4,1e3);
%==============================================%
%% Plots
figure(1)
plot(fr*1e-6, 20*log10(spectr_usrChan),'Linewidth',2)
grid on, grid minor
xlabel("f, MHz")
ylabel("QAM_{out}(f), dB")
set(gca, 'Xlim', [0 400], 'XTick', 0:20:400, 'XTickLabel', 0:20:400)
set(gca, 'fontsize', 15)

figure(2)
plot(t*1e6,sAM,'o-','Linewidth',2)
grid on, grid minor
xlabel("t, us")
ylabel("QAM_{out}(t), V")
set(gca, 'Xlim', [0 1], 'XTick', 0:0.1:1, 'XTickLabel', 0:0.1:1)
set(gca, 'fontsize', 15)

figure(3)
plot(usrDatRsm2I(5:5:end),usrDatRsm2Q(5:5:end),'b.')
hold on
plot(usrDatI(1:5:end)/5,usrDatQ(1:5:end)/5,'or')
grid on, grid minor
xlabel("I")
ylabel("Q")
% set(gca, 'Xlim', [-1 1], 'XTick', -1:0.2:1, 'XTickLabel', -1:0.2:1)
% set(gca, 'Ylim', [-1 1], 'YTick', -1:0.1:1, 'YTickLabel', -1:0.1:1)

figure(4)
plot(usrDatI(1:5:end)/5,usrDatQ(1:5:end)/5,'or')
hold on
plot(usrDatFltI(9:5:end),usrDatFltQ(9:5:end),'b.')

figure(5)
plot(usrDatI(1:5:end)/5,usrDatQ(1:5:end)/5,'or')
hold on
plot(usrDat_Rx_I_symbols,usrDat_Rx_Q_symbols,'b.')
set(gca, 'Xlim', [-1.2 1.2], 'XTick', -1.2:0.2:1.2, 'XTickLabel', -1.2:0.2:1.2)
set(gca, 'Ylim', [-1 0.8], 'YTick', -1:0.2:0.8, 'YTickLabel', -1:0.2:0.8)
grid on, grid minor
xlabel("I")
ylabel("Q")