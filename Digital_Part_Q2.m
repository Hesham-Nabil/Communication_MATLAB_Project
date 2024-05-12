clc;
clear;
close all
bits = 64;
stream =randi([0, 1], 1, bits);

ts = 0.01;
T = bits;
t = 0 : 0.01 : (bits-0.01);
fs = 1 / ts;
df =  1 / T;
f = -0.5 * fs : df : 0.5 * fs - df;
unipolar=zeros(size(t));
for i = 1:bits
    if stream(i) == 1
        unipolar((i-1)*100+1:i*100) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%% Unipolar to freq Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UNIPOLAR= fftshift(fft(unipolar))*ts;
%%%%%%%%%%%%%%%%%%%%%%% Ask Signal Formation in time & Freq Domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=5;
c=3*cos(2*pi*fc*t);
C_30=3*cos(2*pi*fc*t+30);
C_60=3*cos(2*pi*fc*t+60);
C_90=3*cos(2*pi*fc*t+90);
transmitted=c.*unipolar;
Transmitted= fftshift(fft(transmitted))*ts;
%%%%%%%%%%%%%%%%%%%%%%% Plotting Unipolar and Ask Temporal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(2,1,1);
plot(t, unipolar);
xlabel('time');
title('Unipolar in Time Domain');
subplot(2,1,2);
plot(t, transmitted);
xlabel('time');
title('Ask in Time Domain');
%%%%%%%%%%%%%%%%%%%%%%% Plotting Unipolar and Ask Spectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,1,1);
plot(f, abs(UNIPOLAR));
xlabel('Frequency (Hz)');
title('Spectral Domain Of Unipolar');
subplot(2,1,2)
plot(f, abs(Transmitted));
xlabel('Frequency (Hz)');
title('Spectral Domain Of Ask');
grid on;
%%%%%%%%%%%%%%%%%%%%%%% Multiplying Carriers with different phase errors to Ask at Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mixer_output=transmitted.*(c);
mixer_output_30=C_30.*transmitted;
mixer_output_60=C_60.*transmitted;
mixer_output_90=C_90.*transmitted;
%%%%%%%%%%%%%%%%%%%%%%% Comparing in time domain different phase errors After Carrier Multiplication in Mixer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(4,1,1);
plot (t,mixer_output);
xlabel('time');
title('ASK after Carrier Phase 0');
subplot(4,1,2);
plot (t,mixer_output_30);
xlabel('time');
title('ASK after Carrier phase 30');
subplot(4,1,3);
plot (t,mixer_output_60);
xlabel('time');
title('ASK after Carrier phase 60');
subplot(4,1,4);
plot(t, mixer_output_90);
xlabel('time');
title('ASK after Carrier phase 90');
grid on;
%%%%%%%%%%%%%%%%%%%%%%% Comparing in Freq Domain different phase errors After Carrier Multiplication in Mixer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mixer_Output= fftshift(fft(mixer_output))*ts;
Mixer_Output_30= fftshift(fft(mixer_output_30))*ts;
Mixer_Output_60= fftshift(fft(mixer_output_60))*ts;
Mixer_Output_90= fftshift(fft(mixer_output_90))*ts;
figure(4)
subplot(4,1,1);
plot (f,abs(Mixer_Output));
xlabel('frequency');
title('ASK after Carrier phase 0');
subplot(4,1,2);
plot (f,abs(Mixer_Output_30));
xlabel('frequency');
title('ASK after Carrier phase 30');
subplot(4,1,3);
plot (f,abs(Mixer_Output_60));
xlabel('frequency');
title('ASK after Carrier phase 60');
subplot(4,1,4);
plot (f,abs(Mixer_Output_90));
xlabel('frequency');
title('ASK after Carrier phase 90');
%%%%%%%%%%%%%%%%%%%%%%% LPF Applied at f=5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H= abs(f)<= 5;
LPF_Output=H.*Mixer_Output;
LPF_Output_30=H.*Mixer_Output_30;
LPF_Output_60=H.*Mixer_Output_60;
LPF_Output_90=H.*Mixer_Output_90;
%%%%%%%%%%%%%%%%%%%%%%% IFT to get Filtered Ask Received Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
received =real(ifft(ifftshift(LPF_Output)/ts));
received_30 =real(ifft(ifftshift(LPF_Output_30)/ts));
received_60 =real(ifft(ifftshift(LPF_Output_60)/ts));
received_90 =real(ifft(ifftshift(LPF_Output_90)/ts));
figure(5)
subplot(4,1,1);
plot (t,received);
xlabel('time');
title('ASK phase 0 Received'); %Received Ask Message Without Phase Error
subplot(4,1,2);
plot (t,received_30);
xlabel('time');
title('ASK phase 30 Received'); %Phase 30 Amplitudes are Aprroximately Equal to Received Ask Amplitudes divided by 6.4 , Also Lowest SNR
subplot(4,1,3);
plot (t,received_60);
xlabel('time');
title('ASK phase 60 Received'); %Phase 60 Amplitudes are Equal to Received Ask Amplitudes Multiplied by -1 (Inverted)
subplot(4,1,4);
plot (t,received_90);
xlabel('time');
title('ASK phase 90 Received');  %Phase 90 Amplitudes are Approximately Equal to Phase 60 Amplitudes divided by 2
output= zeros(size(t));
for i = 1:(bits*100)
    if received(i) <0.5 && received(i) >-0.5 
        output(i) = 0;
    else
        output(i) = 1;
    end
end
figure(6)
plot (t,output);
xlabel('time');
title('ASK Received');
output_30= zeros(size(t));
for i = 1:(bits*100)
    if received_30(i) <0.5 && received_30(i) >-0.5 
        output_30(i) = 0;
    else
        output_30(i) = 1;
    end
end
figure(7)
plot (t,output_30);
xlabel('time');
title('ASK phase 30 Received');

output_60= zeros(size(t));
for i = 1:(bits*100)
    if received_60(i) <0.5 && received_60(i) >-0.5 
        output_60(i) = 0;
    else
        output_60(i) = 1;
    end
end
figure(8)
plot (t,output_60);
xlabel('time');
title('ASK phase 60 Received');

output_90= zeros(size(t));
for i = 1:(bits*100)
    if received_90(i) <0.5 && received_90(i) >-0.5 
        output_90(i) = 0;
    else
        output_90(i) = 1;
    end
end
figure(9)
plot (t,output_90);
xlabel('time');
title('ASK phase 600 Received');
