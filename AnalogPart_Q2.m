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
m(t >= 0 & t < 6) = cos( 2*pi*t(t >= 0 & t < 6));     % m(t) = cos(2*pi*t) 0<= t <= 6


figure(1)
plot(t,m);
xlabel('time');
ylabel('Baseband Signal');
title('Plot Of X(t)');

#################  frequency respone of base band signal  #################

M = fftshift(fft(m))*ts;
Analytical = 3*sinc(6*(f-1))+3*sinc(6*(f+1));

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
title('3*sinc(6*f-1)+3*sinc(6*f+1)');


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
