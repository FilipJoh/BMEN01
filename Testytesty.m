%% load data and define constants
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;

%% Decimate signal
df = 2;
ecg_h = decimate(ecg_h,df);
t_h = round(t_h / df);
F_s = F_s / df;
T = 1 / F_s;

%% Filter signal(s)
% cutFreq = 50; %Hz
% B = fir1(10000,cutFreq/(F_s/2));
% [H,w]=freqz(B,1);
% plot(F_s * w/2/pi,abs(H))
% ecg_h = filter(B,1,ecg_h);

%% Pair consecutive heartbeats
pairs = 3;
offset = 0;
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
    
    duration_per_heartbeat = T * L / 2;
    xval = (1 / duration_per_heartbeat)/2;
    x=[xval,xval];
    y=[0,1000];
    plot(x,y,'-r');
    plot(-x,y,'-r');
end

%% TEST
nbrOfWindows = 1;
offset = 500;
figure;
for i=1:(nbrOfWindows)
    k = 2*i - 1;
    signal = [];
    nbrOfBeats = 2;
    for a = 1:nbrOfBeats
        signal = [signal ecg_h((t_h(k + a-1 + offset):t_h(k + a-1 + offset + 1)))];
    end
    L = length(signal);
    
    duration_per_heartbeat = T * L / nbrOfBeats;
    xval = (1 / duration_per_heartbeat)/2;
    
    
    subplot(2,nbrOfWindows,i);
    plot(T * (1:L),signal);
    xlabel('time /s');
    ylabel('amplitude /mV');
    axis([0 T*L -inf inf]);
    
    NFFT = 2^(nextpow2(L) + 5);
%     w = hanning(L,'periodic');
%     y = signal;
%     y = w'.*signal;
%     LSD = abs(fft(y,NFFT)) / L;
%     f = F_s * linspace(0,1,NFFT);
    [LSD,f] = pyulear(signal,L,NFFT,F_s);
%     [Y,f] = periodogram(signal,rectwin(L), NFFT,F_s);
    
    subplot(2,nbrOfWindows,i + nbrOfWindows);
    LSD_log = 10*log10(LSD);
    plot(f,LSD,'.');
    xlabel('Frequency /Hz');
    ylabel('PSD / dB');
    axis([0 5 0 1.2 * max(LSD)]);
    
    xmark=[xval,xval];
    ymark=[0,10^10];
    hold on;
    plot(xmark,ymark,'-r');
end