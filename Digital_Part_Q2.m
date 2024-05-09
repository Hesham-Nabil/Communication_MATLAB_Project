clc;
clear;
close all
bits = 64;
stream =randi([0, 1], 1, bits);

ts = 0.05;
T = bits;
t = 0 : 0.01 : (bits-0.01);
fs = 1 / ts;
df =  0.01*fs / T;
f = -0.5 * fs : df : 0.5 * fs - df;
unipolar=zeros(size(t));
for i = 1:bits
    if stream(i) == 1
        unipolar((i-1)*100+1:i*100) = 1;
    end
end

UNIPOLAR= fftshift(fft(unipolar))*ts;
fc=5;
c=3*cos(2*pi*fc*t);
C_30=3*cos(2*pi*fc*t+30);
C_60=3*cos(2*pi*fc*t+60);
C_90=3*cos(2*pi*fc*t+90);
Ask=c.*unipolar;
Ask_frequency= fftshift(fft(Ask))*ts;
figure(1);
subplot(2,1,1);
plot(t, unipolar);
xlabel('time');
title('Unipolar in Time Domain');
subplot(2,1,2);
plot(t, Ask);
xlabel('time');
title('Ask in Time Domain');
figure(2)
subplot(2,1,1);
plot(f, abs(UNIPOLAR));
xlabel('Frequency (Hz)');
title('Spectral Domain Of Unipolar');
subplot(2,1,2)
plot(f, abs(Ask_frequency));
xlabel('Frequency (Hz)');
title('Spectral Domain Of Ask');
grid on;
Mixer_Output_Time=Ask.*(c);
Mixer_Output_Time_30=C_30.*Ask;
Mixer_Output_Time_60=C_60.*Ask;
Mixer_Output_Time_90=C_90.*Ask;
figure(3)
subplot(4,1,1);
plot (t,Mixer_Output_Time);
xlabel('time');
title('ASK after Carrier Phase 0');
subplot(4,1,2);
plot (t,Mixer_Output_Time_30);
xlabel('time');
title('ASK after Carrier phase 30');
subplot(4,1,3);
plot (t,Mixer_Output_Time_60);
xlabel('time');
title('ASK after Carrier phase 60');
subplot(4,1,4);
plot(t, Mixer_Output_Time_90);
xlabel('time');
title('ASK after Carrier phase 90');
grid on;
Mixer_Output_Freq= fftshift(fft(Mixer_Output_Time))*ts;
Mixer_Output_Freq_30= fftshift(fft(Mixer_Output_Time_30))*ts;
Mixer_Output_Freq_60= fftshift(fft(Mixer_Output_Time_60))*ts;
Mixer_Output_Freq_90= fftshift(fft(Mixer_Output_Time_90))*ts;
figure(4)
subplot(4,1,1);
plot (f,abs(Mixer_Output_Freq));
xlabel('frequency');
title('ASK after Carrier phase 0');
subplot(4,1,2);
plot (f,abs(Mixer_Output_Freq_30));
xlabel('frequency');
title('ASK after Carrier phase 30');
subplot(4,1,3);
plot (f,abs(Mixer_Output_Freq_60));
xlabel('frequency');
title('ASK after Carrier phase 60');
subplot(4,1,4);
plot (f,abs(Mixer_Output_Freq_90));
xlabel('frequency');
title('ASK after Carrier phase 90');
H= abs(f)<= 1;
LPF_Output=H.*Mixer_Output_Freq;
LPF_Output_30=H.*Mixer_Output_Freq_30;
LPF_Output_60=H.*Mixer_Output_Freq_60;
LPF_Output_90=H.*Mixer_Output_Freq_90;
Ask_Received =real(ifft(ifftshift(LPF_Output)/ts));
Ask_Received_30 =real(ifft(ifftshift(LPF_Output_30)/ts));
Ask_Received_60 =real(ifft(ifftshift(LPF_Output_60)/ts));
Ask_Received_90 =real(ifft(ifftshift(LPF_Output_90)/ts));
figure(5)
subplot(4,1,1);
plot (t,Ask_Received);
xlabel('time');
title('ASK phase 0 Received');
subplot(4,1,2);
plot (t,Ask_Received_30);
xlabel('time');
title('ASK phase 30 Received');
subplot(4,1,3);
plot (t,Ask_Received_60);
xlabel('time');
title('ASK phase 60 Received');
subplot(4,1,4);
plot (t,Ask_Received_90);
xlabel('time');
title('ASK phase 90 Received');
Output=zeros(size(t));
for i = 1:6400
    if (Ask_Received(i)>=2.25) Output(i)=1;
    end
end
Output_30=zeros(size(t));
for i = 1:6400
    if (Ask_Received_30(i)>=2.25) Output_30(i)=1;
    end
end
Output_60=zeros(size(t));
for i = 1:6400
    if (Ask_Received_60(i)>=2.25) Output_60(i)=1;
    end
end
Output_90=zeros(size(t));
for i = 1:6400
    if (Ask_Received_90(i)>=2.25) Output_90(i)=1;
    end
end
figure(6)
subplot(4,1,1);
plot (t,Output);
xlabel('time');
title('Output Phase 0');
subplot(4,1,2);
plot (t,Output_30);
xlabel('time');
title('Output Phase 30');
subplot(4,1,3);
plot (t,Output_60);
xlabel('time');
title('Output Phase 60');
subplot(4,1,4);
plot (t,Output_90);
xlabel('time');
title('Output Phase 90');