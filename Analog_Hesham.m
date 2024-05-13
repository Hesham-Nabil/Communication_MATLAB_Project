clc;
clear;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%    Time Domain    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ts = 0.01;                                          % Time step
T = 100;                                            % simulation time
N = ceil(T/ts);                                     % number of time samples
t = -(N/2)*ts : ts : ((N-1)/2)*ts;                  % time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%  frequency domain  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 100;                                           % 1/ts
df = 0.01;                                          % fs/N
if (rem(N,2)==0)
  f = -(0.5*fs): df : (0.5*fs-df);
else
  f = -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  BaseBand Signal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(size(t));
x(t >= -2 & t < -1) = t(t >= -2 & t < -1)+2;        % m(t) = t+2 for the range -2 <= t <=-1
x(t >= -1 & t <= 1) = 1;                            % m(t) = 1   for the range -1 <= t <= 1
x(t > 1 & t <= 2) = -t(t > 1 & t <= 2)+2;           % m(t) =-t+2 for the range  1 <= t <= 2

figure(1)                                           % Plot Of X(t)
plot(t,x);
xlabel('time');
ylabel('Baseband Signal');
title('Plot Of X(t)');

%%%%%%%%%%%%%%%%%  frequency respone of base band signal  %%%%%%%%%%%%%%%%%

X = fftshift(fft(x))*ts;                             % Practical Expression
Analytical = -3 * sinc(f) .* sinc(3*f);              % Analytical Expression

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(1,2,1)
plot(f,abs(X));
xlabel('time');
ylabel('Frequency respone of x(t)');
title('X(f)');
subplot(1,2,2)
plot(f,abs(Analytical));
xlabel('frequency');
ylabel('Frequency respone of x(t)');
title('Analytical Solution');


%%%%%%%%%%%%%%%%%%%%%%%%        Bandwidth       %%%%%%%%%%%%%%%%%%%%%%%%%%%

Total_Energy_in_Freq = sum(abs(X).^2)*df;      % ∑ M^2 df
zero_freq = find(f==0);                        % since the spectrum is even around the y-axis, 
                                               % we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(f) )
  Energy_accumulator =  Energy_accumulator + (abs(X(index)).^2)*df;     % accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq)               % 0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    Bandwidth = f(index)
    break
  end
end

%%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=1HZ   %%%%%%%%%%%%%%%%%%%%

H1 = zeros(size(f));
H1 = abs(f)<(1+df);
Filtered_Signal1 = X.*H1;
filtered_signal1= real(ifft(ifftshift(Filtered_Signal1) /ts));

figure(3)
subplot(2,2,1)
plot(f,abs(H1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 1 Hz')
subplot(2,2,2)
plot(f,abs(Filtered_Signal1))
xlabel('Frequency (Hz)')
ylabel('Signal')
title('Filtered_Signal with BW = 1Hz')
subplot(2,1,2)
plot(t,filtered_signal1,'r',t,x,'b');
xlabel('Time')
ylabel('Signal')
title('Filtered_Signal with BW = 1Hz')
legend('Input Signal','Filtered Signal')

%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=0.3HZ   %%%%%%%%%%%%%%%%%%%

H2 = zeros(size(Analytical));
H2 = abs(f)<(0.3+df);
Filtered_Signal2 = X.*H2;
filtered_signal2= real(ifft(ifftshift(Filtered_Signal1) /ts));

figure(4)
subplot(2,2,1)
plot(f,abs(H2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 0.3 Hz')
subplot(2,2,2)
plot(f,abs(Filtered_Signal2))
xlabel('Frequency (Hz)')
ylabel('Signal')
title('Filtered_Signal with BW = 0.3 Hz')
subplot(2,1,2)
plot(t,filtered_signal2,'r',t,x,'b');
xlabel('Time')
ylabel('Signal')
title('Filtered_Signal with BW = 0.3Hz')
legend('Input Signal','Filtered Signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  BaseBand Signal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = zeros(size(t));
m(t >= 0 & t< 6) = cos( 2*pi*t(t >= 0 & t < 6));  % m(t) = cos(2*pi*t) 0<= t <= 6
figure(7)
plot(t,m);
xlabel('time(sec)');
ylabel('Baseband Signal m(t)');
title('Plot Of m(t)');

%%%%%%%%%%%%%%%%%  frequency respone of base band signal  %%%%%%%%%%%%%%%%%

M = fftshift(fft(m))*ts;
Analyticalm = 3*sinc(6*(f-1))+3*sinc(6*(f+1)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
subplot(1,2,1)
plot(f,abs(M));
xlabel('Frequency Hz');
ylabel('Frequency respone of m(t)');
title('M(f)');
subplot(1,2,2)
plot(f,abs(Analyticalm));
xlabel('Frequency HZ');
ylabel('Frequency respone of m(t)');
title('Analytical Solution');

%%%%%%%%%%%%%%%%%%%%%%%%%        Bandwidth       %%%%%%%%%%%%%%%%%%%%%%%%%%

Total_Energy_in_Freq = sum(abs(M).^2)*df;                  % ∑ M^2 df

zero_freq = find(f==0);                                                     
Energy_accumulator=0;
for(index = zero_freq : length(f) )
  Energy_accumulator =  Energy_accumulator + (abs(M(index)).^2)*df;  
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); 
    Bandwidthm = f(index)
    break
  end
end

%%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=1HZ   %%%%%%%%%%%%%%%%%%%%

H1 = zeros(size(f));
H1 = abs(f)<(1+df);
Filtered_Signal3 = M.*H1;
filtered_signal3= real(ifft(ifftshift(Filtered_Signal3) /ts));

figure(3)
subplot(2,2,1)
plot(f,abs(H1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 1 Hz')
subplot(2,2,2)
plot(f,abs(Filtered_Signal3))
xlabel('Frequency (Hz)')
ylabel('Signal')
title('Filtered_Signal with BW = 1Hz')
subplot(2,1,2)
plot(t,filtered_signal3,'r',t,m,'b');
xlabel('Time')
ylabel('Signal')
title('Filtered_Signal with BW = 1Hz')
legend('Input Signal','Filtered Signal')

%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=0.3HZ   %%%%%%%%%%%%%%%%%%%

H2 = zeros(size(Analytical));
H2 = abs(f)<(0.3+df);
Filtered_Signal4 = M.*H2;
filtered_signal4= real(ifft(ifftshift(Filtered_Signal4) /ts));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
subplot(2,2,1)
plot(f,abs(H2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 0.3 Hz')
subplot(2,2,2)
plot(f,abs(Filtered_Signal4))
xlabel('Frequency (Hz)')
ylabel('Signal')
title('Filtered_Signal with BW = 0.3 Hz')
subplot(2,1,2)
plot(t,filtered_signal4,'r',t,m,'b');
xlabel('Time')
ylabel('Signal')
title('Filtered_Signal with BW = 0.3Hz')
legend('Input Signal','Filtered Signal')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% FDM MODULATION SCHEME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Signal 1 (DSB-SC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc1=20;                                                 % Carrier Freq 1
c1= cos(2*pi*fc1*t);                                    % Carrier Signal 1
s1 = filtered_signal1 .*c1 ;                             % modulated x(t)
S1= fftshift(fft(s1))*ts;                               % Transform of S1
F1= zeros(size(s1));                                    % BPF FOR X(T)
F1= ((abs(f) > (18.5-df )) & ((21.5+df) > abs(f)));     % FROM 18.5 TO 21.5 FOR X(T)
tx_s1 = real(ifft(ifftshift(S1.*F1) / ts));    % s1 in time domain

%%%%%%%%%%%%%%%%%%%%%%%%%%   Signal 2 (SSB) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc2=25.5;                                               % Carrier Freq 2
c2= cos(2*pi*fc2*t);                                    % Carrier Signal 2

s2 = m .*c2;                                            % modulated m(t)
S2=fftshift(fft(s2))*ts;                                % Transform of S2

F2= zeros(size(s2));                                    % BPF FOR M(T)
F2 = ((abs(f) > (25.5-df)) & ((27 + df) > abs(f)));     % FROM 25.5 TO 27 FOR SSB

tx_s2 = real(ifft(ifftshift(S2.*F2) / ts));    % s2 in time domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% S1 + S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx_s=tx_s1+tx_s2;            %  s1 + S2 in time domain
TX_S= fftshift(fft(tx_s))*ts;         %  frequency spectrum of (s1 + s2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)       
subplot(2,1,1)
plot(t,tx_s);
xlabel('Time (sec)');
ylabel('Time respone of Transmitted s(t)');
title(' s1(t) + s2(t) ');
subplot(2,1,2)
plot(f,abs(TX_S));
xlabel('Frequency HZ');
ylabel('Frequency respone of Transmitted s(t)');
title(' s1(F) + s2(F) ');

%%%%%%%%%%%%%%%%%%%%%%%%  Coherent Demodulator  %%%%%%%%%%%%%%%%%%%%%%%%%%%

Filter= zeros(size(TX_S));                      % LPF for demodulating
Filter= abs(f)<=1.5+df;                                  % 1.5 as band width = 3

RX_S1=fftshift(fft(tx_s.*c1))*ts.*Filter; % recieved S1 AFTER LPF
rx_s1 = real(ifft(ifftshift(RX_S1) /ts));    % s1 in time domain

RX_S2=fftshift(fft(tx_s.*c2))*ts.*Filter; % Recieved s2 after lpf
rx_s2 = real(ifft(ifftshift(RX_S2) /ts));    % s2 in time domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
subplot(2,2,1)
plot(t, rx_s1/max(rx_s1), 'b-', t, filtered_signal1/max(filtered_signal1), 'r-');
legend("Recieved Signal","Input Signal");
xlabel('Time (sec)');
ylabel('Normalized Time respone of Recieved and input signals');
title(' s1(t)/filtered signal(t) ');
subplot(2,2,2)
plot(t, rx_s2/max(rx_s2), 'b-', t, m/max(m), 'r-');
legend("Recieved Signal","Input Signal");
xlabel('Time (sec)');
ylabel('Normalized Time respone of Recieved and input signals');
title(' s2(t)/m(t) ');
subplot(2,2,3)
plot(f, abs(RX_S1)/max(abs(RX_S1)), 'b-', f, abs(Filtered_Signal1)/max(abs(Filtered_Signal1)), 'r-');
legend("Recieved Signal","Input Signal");
xlabel('Frequency (HZ)');
ylabel('Normalized Frequency respone of Recieved and input signals');
title(' S1(F)/Filtered Signal (F) ');
subplot(2,2,4)
plot(f, abs(RX_S2)/max(abs(RX_S2)), 'b-', f, abs(M)/max(abs(M)), 'r-');
legend("Recieved Signal","Input Signal");
xlabel('Frequency (HZ)');
ylabel('Normalized Frequency respone of Recieved and input signals');
title(' S2(F)/M(F) ');
