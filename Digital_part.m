clc;
clear;
close all
bits = 70;
% stream =randi([0, 1], 1, bits);S
% stream = zeros(1,70);
stream = ones(1,70);
% stream(1)=1;
ts = 0.01;
T = bits;
t = 0 : 0.01 : (bits-0.01);
fs = 1 / ts;
df =  1 / T;
f = -0.5 * fs : df : 0.5 * fs - df;

bipolar = zeros(size(t));
flag = 1;  % Start with +ve pulse
for i = 1:bits
    if stream(i) == 1
       flag=-flag ;
       bipolar((i-1)*100+1:i*100) = flag;
    end
end

unipolar=zeros(size(t));
for i = 1:bits
    if stream(i) == 1
        unipolar((i-1)*100+1:i*100) = 1;
    end
end
BIPOLAR = fftshift(fft(bipolar))*ts;
UNIPOLAR= fftshift(fft(unipolar))*ts;

Total_Energy_in_Freq = sum(abs(BIPOLAR).^2)*df;      % ∑ M^2 df
zero_freq = find(f==0);   %since the spectrum is even around the y-axis, we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(f) )
  Energy_accumulator =  Energy_accumulator + (abs(BIPOLAR(index)).^2)*df;   %accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); %0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    BandwidthB = f(index)
    break
  end
end
Total_Energy_in_Freq = sum(abs(UNIPOLAR).^2)*df; 
zero_freq = find(f==0);   %since the spectrum is even around the y-axis, we need to find the component number in the middle at f==0
Energy_accumulator=0;
for(index = zero_freq : length(f) )
  Energy_accumulator =  Energy_accumulator + (abs(UNIPOLAR(index)).^2)*df;   %accumulating the energy of all components until we reach the end of the spectrum
  if(Energy_accumulator >= (0.95/2)*Total_Energy_in_Freq); %0.95/2 because we are calculating the power of the right half of the spectrum so the total power is expected to be 50%
    BandwidthU = f(index)
    break
  end
end

figure(1)
plot(t, bipolar);
xlabel('Time');
ylabel('Amplitude');
title('Bipolar');
figure(2)
plot(t, unipolar);
xlabel('Time');
ylabel('Amplitude');
title('Unipolar');
figure(3)
plot(f, abs(BIPOLAR));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectral Domain Of Bipolar');
figure(4)
plot(f, abs(UNIPOLAR));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain Of Unipolar');
grid on;
