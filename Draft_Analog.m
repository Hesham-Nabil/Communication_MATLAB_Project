clc;
clear;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%    Time Domain    %%%%%%%%%%%%%%%%%%%%%%%%%%%

tsx = 0.01;                   %Time step
Tx = 100;                     %simulation time
Nx = ceil(Tx/tsx);              %number of time samples
tx = -(Nx/2)*tsx : tsx : ((Nx-1)/2)*tsx;       %time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%  frequency domain  %%%%%%%%%%%%%%%%%%%%%%%%%%

fsx = 100;   %1/ts
dfx = 0.01;  %fs/N
if (rem(Nx,2)==0)
  fx = -(0.5*fsx): dfx : (0.5*fsx-dfx);
else
  fx = -(0.5*fsx-0.5*dfx): dfx : (0.5*fsx-0.5*dfx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  BaseBand Signal  %%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(size(tx));
x(tx >= -2 & tx < -1) = tx(tx >= -2 & tx < -1)+2;     % m(t) = t+2 for the range -2 <= t <= -1
x(tx >= -1 & tx <= 1) = 1;                         % m(t) = 1   for the range -1 <= t <= 1
x(tx > 1 & tx <= 2) = -tx(tx > 1 & tx <= 2)+2;        % m(t) =-t+2 for the range  1 <= t <= 2

figure(1)
plot(tx,x);
xlabel('time');
ylabel('Baseband Signal');
title('Plot Of X(t)');

%%%%%%%%%%%%%%%%%  frequency respone of base band signal  %%%%%%%%%%%%%%%%%

X = fftshift(fft(x))*tsx;
Analyticalx = -3 * sinc(fx) .* sinc(3*fx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
subplot(1,2,1)
plot(fx,abs(X));
xlabel('time');
ylabel('Frequency respone of x(t)');
title('X(f)');
subplot(1,2,2)
plot(fx,abs(Analyticalx));
xlabel('frequency');
ylabel('Frequency respone of x(t)');
title('Analytical Solution');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%        Bandwidth       %%%%%%%%%%%%%%%%%%%%%%%%%

Total_Energy_in_Freq = sum(abs(X).^2)*dfx;      % ∑ M^2 df

zero_freq = find(fx==0);   %since the spectrum is even around the y-axis, we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(fx) )
  Energy_accumulator =  Energy_accumulator + (abs(X(index)).^2)*dfx;   %accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); %0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    Bandwidthx = fx(index)
    break
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=1HZ   %%%%%%%%%%%%%%%%%%%%%%%
H1 = zeros(size(Analyticalx));
H1 = abs(fx)<(1+dfx);
Filtered_Signal1 = X.*H1;

figure(2)
subplot(1,2,1)
plot(fx,abs(H1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 1 Hz')
subplot(1,2,2)
plot(fx,abs(Filtered_Signal1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filtered_Signal with BW = 1Hz')

%%%%%%%%%%%%%%%%%%%%%%%%%   Ideal LowPass Filter BW=0.3HZ   %%%%%%%%%%%%%%%%%%%%%%

H2 = zeros(size(Analyticalx));
H2 = abs(fx)<(0.3+dfx);
Filtered_Signal2 = X.*H2;

figure(3)
subplot(1,2,1)
plot(fx,abs(H2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 0.3 Hz')
subplot(1,2,2)
plot(fx,abs(Filtered_Signal2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filtered_Signal with BW = 0.3 Hz')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%    Time Domain    %%%%%%%%%%%%%%%%%%%%%%%%%%%

tsm = 0.01;                   %Time step
Tm = 100;                     %simulation time
Nm = ceil(Tm/tsm);              %number of time samples
tm = -(Nm/2)*tsm : tsm : ((Nm-1)/2)*tsm;       %time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%  frequency domain  %%%%%%%%%%%%%%%%%%%%%%%%%%

fsm = 100;   %1/ts
dfm = 0.01;  %fs/N
if (rem(Nm,2)==0)
  fm = -(0.5*fsm): dfm : (0.5*fsm-dfm);
else
  fm = -(0.5*fsm-0.5*dfm): dfm : (0.5*fsm-0.5*dfm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  BaseBand Signal  %%%%%%%%%%%%%%%%%%%%%%%%%%%

m = zeros(size(tm));
m(tm >= 0 & tm< 6) = cos( 2*pi*tm(tm >= 0 & tm < 6));     % m(t) = cos(2*pi*t) 0<= t <= 6


figure(4)
plot(tm,m);
xlabel('time');
ylabel('Baseband Signal');
title('Plot Of m(t)');

%%%%%%%%%%%%%%%%%  frequency respone of base band signal  %%%%%%%%%%%%%%%%%

M = fftshift(fft(m))*tsm;
Analyticalm = 3*sinc(6*(fm-1))+3*sinc(6*(fm+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
subplot(1,2,1)
plot(fm,abs(M));
xlabel('time');
ylabel('Frequency respone of m(t)');
title('M(f)');
subplot(1,2,2)
plot(fm,abs(Analyticalm));
xlabel('frequency');
ylabel('Frequency respone of m(t)');
title('3*sinc(6*f-1)+3*sinc(6*f+1)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%        Bandwidth       %%%%%%%%%%%%%%%%%%%%%%%%%

Total_Energy_in_Freq = sum(abs(M).^2)*dfm;      % ∑ M^2 df

zero_freq = find(fm==0);   %since the spectrum is even around the y-axis, we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(fm) )
  Energy_accumulator =  Energy_accumulator + (abs(M(index)).^2)*dfm;   %accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); %0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    Bandwidthm = fm(index)
    break
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FDM MODULATION SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc1=20;
fc2=25.5;
c1= cos(2*pi*fc1*tx);
c2=cos(2*pi*fc2*tm);
x = real(ifft(ifftshift(Filtered_Signal1) /tsx));%%signal x(t)
s1 = x .*c1;%% modulated x(t)
s2 = m .*c2;%%modulated m(t)
F1= zeros(size(s1));%% BPF FOR M(T)
F2= zeros(size(s2));%%BPF FOR M(T)
F1 = ((abs(fx) > (18.5 )) & ((21.5) > abs(fx)));%FROM 18.5 TO 21.5 FOR X(T)
F2 = ((abs(fm) > (25.5 - dfm)) & ((27 + dfm) > abs(fm)));%FROM 25.5 TO 27 FOR SSB
Transmitted_s1 = real(ifft(ifftshift(fftshift(fft(s1))*tsx.*F1) / tsx));%%s1 in time domain
Transmitted_s2 = real(ifft(ifftshift(fftshift(fft(s2))*tsm.*F2) / tsm));%% s2 in time domain
Transmitted_s=Transmitted_s1+Transmitted_s2;% s in the medieum
Transmitted_S= fftshift(fft(Transmitted_s))*tsx;%frequency spectrum of S
figure(6)
subplot(1,2,1)
plot(tx,Transmitted_s);
xlabel('Time');
ylabel('Time respone of Transmitted_s(t)');
title(' s1(t) + s2(t) ');
subplot(1,2,2)
plot(fx,abs(Transmitted_S));
xlabel('frequency');
ylabel('Frequency respone of Transmitted_s(t)');
title(' s1(F) + s2(F) ');
Filter= zeros(size(Transmitted_S));% LPF for demodulating
Filter= abs(fm)<=1.5+dfx;% 1.5 as band width = 3
Recieved_S1=fftshift(fft(Transmitted_s.*c1))*tsx.*Filter; %%recieved S1 AFTER LPF
Recieved_S2=fftshift(fft(Transmitted_s.*c2))*tsm.*Filter;%%% Recieved s2 after lpf
Recieved_s1 = real(ifft(ifftshift(Recieved_S1) /tsx));  %s1 in time domain
Recieved_s2 = real(ifft(ifftshift(Recieved_S2) /tsm));  % s2 in time domain
figure(7)
plot(tm, Recieved_s2, 'b-', tm, m, 'r-');
%%Recieved_S2=fftshift(fft(Transmitted_s.*c2))*tsm.*LPF;
%%Recieved_s1 = real(ifft(ifftshift(Recieved_S1) * Nx));
%%Recieved_s2 = real(ifft(ifftshift(Recieved_S2) * Nm));
