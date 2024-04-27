clc;
clear;
close all
bits = 64;
stream =randi([0, 1], 1, bits);

ts = bits;
T = 100*bits;
t = 0 : 1 : (100*bits-1);
fs = 1 / ts;
df =  fs / T;
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

figure(1)
subplot(2, 1, 1);
plot(t, bipolar);
xlabel('Time');
ylabel('Amplitude');
title('Time Domain');
subplot(2, 1, 2);
plot(t, unipolar);
xlabel('Time');
ylabel('Amplitude');
title('Unipolar');
BIPOLAR = fftshift(fft(bipolar))*ts;
UNIPOLAR= fftshift(fft(unipolar))*ts;
figure(2)
subplot(2,1,1);
plot(f, abs(BIPOLAR));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectral Domain Of Bipolar');
subplot(2,1,2)
plot(f, abs(UNIPOLAR));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectral Domain Of Unipolar');
grid on;
