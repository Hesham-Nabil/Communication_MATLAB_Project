clc;
clear;
close all
bits = 64;
stream =randi([0, 1], 1, bits);
% stream = zeros(1,64);
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
BIPOLAR = fftshift(fft(bipolar))*ts;
UNIPOLAR= fftshift(fft(unipolar))*ts;

power_spectrum = abs(BIPOLAR).^2;  % Power spectrum

% Determine bandwidth
max_power = max(power_spectrum);  % Maximum power
threshold_power = 0.05 * max_power;  % Threshold power (5% of max power)
bandwidth_indices = find(power_spectrum >= threshold_power);  % Find indices where power exceeds threshold
bandwidth = abs(f(bandwidth_indices(end)) - f(bandwidth_indices(1)));  % Compute bandwidth

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
