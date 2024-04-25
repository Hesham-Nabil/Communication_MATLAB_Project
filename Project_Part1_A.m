clc;
clear;
close all

###########################    Time Domain    ###########################

ts = 0.01;                   %Time step
T = 100;                     %simulation time
N = ceil(T/ts);              %number of time samples
t = -(N/2)*ts : ts : ((N-1)/2)*ts;       %time vector

###########################  frequency domain  ##########################

fs = 100;   %1/ts
df = 0.01;  %fs/N
if (rem(N,2)==0)
  f = -(0.5*fs): df : (0.5*fs-df);
else
  f = -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df);
end

###########################  BaseBand Signal  ###########################

m = zeros(size(t));
m(t >= -2 & t < -1) = t(t >= -2 & t < -1)+2;     % m(t) = t+2 for the range -2 <= t <= -1
m(t >= -1 & t <= 1) = 1;                         % m(t) = 1   for the range -1 <= t <= 1
m(t > 1 & t <= 2) = -t(t > 1 & t <= 2)+2;        % m(t) =-t+2 for the range  1 <= t <= 2

figure(1)
plot(t,m);
xlabel('time');
ylabel('Baseband Signal');
title('Plot Of X(t)');

#################  frequency respone of base band signal  #################

M = fftshift(fft(m))*ts;
Analytical = -3 * sinc(f) .* sinc(3*f);

##############################    plotting   ##############################

figure(2)
subplot(1,2,1)
plot(f,abs(M));
xlabel('time');
ylabel('Frequency respone of m(t)');
title('M(f)');
subplot(1,2,2)
plot(f,abs(Analytical));
xlabel('frequency');
ylabel('Frequency respone of m(t)');
title('Analytical Solution');


############################        Bandwidth       #########################

Total_Energy_in_Freq = sum(abs(M).^2)*df;      % ∑ M^2 df

zero_freq = find(f==0);   %since the spectrum is even around the y-axis, we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(f) )
  Energy_accumulator =  Energy_accumulator + (abs(M(index)).^2)*df;   %accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); %0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    Bandwidth = f(index)
    break
  end
end

##########################   Ideal LowPass Filter BW=1HZ   #######################
H1 = zeros(size(Analytical));
H1 = abs(f)<(1+df);
Filtered_Signal1 = M.*H1;

figure(2)
subplot(1,2,1)
plot(f,abs(H1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 1 Hz')
subplot(1,2,2)
plot(f,abs(Filtered_Signal1))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filtered_Signal with BW = 1Hz')

#########################   Ideal LowPass Filter BW=0.3HZ   ######################

H2 = zeros(size(Analytical));
H2 = abs(f)<(0.3+df);
Filtered_Signal2 = M.*H2;

figure(3)
subplot(1,2,1)
plot(f,abs(H2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filter Spectrum With BW = 0.3 Hz')
subplot(1,2,2)
plot(f,abs(Filtered_Signal2))
xlabel('Frequency (Hz)')
ylabel('|H(f)|')
title('Filtered_Signal with BW = 0.3 Hz')



