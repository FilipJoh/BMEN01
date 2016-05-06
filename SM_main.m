%% load data and define constants
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;

durationsSamples = (t_h(2:end) - t_h(1:(end-1)));
durations = T * durationsSamples;
meanDuration = mean(durations);
meanHeartrate = 60*(1 / meanDuration);

%% Filter signal(s)
cutFreq = 50; %Hz
B = fir1(10000,cutFreq/(F_s/2));
[H,w]=freqz(B,1);
plot(F_s * w/2/pi,abs(H))
ecg_h = filter(B,1,ecg_h);

%% Visualize heartbeat pairs vertically
pairs = 3;
offset = 599;
figure;
for i=1:(pairs)
    k = 2*i - 1;
    range1 = (t_h(k + offset):t_h(k + offset + 1));
    range2 = (t_h(k + 1 + offset):t_h(k + 1 + offset + 1));
    heartbeat1 = ecg_h(range1);
    heartbeat2 = ecg_h(range2);
    
    signal = [heartbeat1 heartbeat2];
    L = length(signal);
    
    subplot(2,pairs,i);
    plot(T * (1:L),signal);
    xlabel('time /s');
    ylabel('amplitude /mV');
    axis([0 T*L -inf inf]);
    
    NFFT = 2^nextpow2(L)*2^10;
    Y = fft(signal,NFFT) / L;
    Y = abs(fftshift(Y));
    f = F_s * linspace(-0.5,0.5,NFFT);
    
    subplot(2,pairs,i + pairs);
    plot(f,Y)
    xlabel('Frequency /Hz')
    ylabel('|Y(f)|')
    axis([0 5 0 1.2 * max(Y)]);
    hold on;
    
    duration_local = T * L / 2;
    xval = (1 / duration_local)/2;
    x=[xval,xval];
    y=[0,1000];
    plot(x,y,'-r');
    plot(-x,y,'-r');
end