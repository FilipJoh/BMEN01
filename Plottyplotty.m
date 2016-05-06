%% load data and define constants
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;
t = T * 1:length(ecg_h);

% plot(t,ecg_h);

durationsSamples = (t_h(2:end) - t_h(1:(end-1)));
durations = T * durationsSamples;
meanDuration = mean(durations);
meanHeartrate = 60*(1 / meanDuration);

%% Filter signal(s)
% cutFreq = 60; %Hz
% B = fir1(100,cutFreq/(F_s/2));
% [H,w]=freqz(B,1);
% plot(F_s * w/2/pi,abs(H))
% ecg_h = filter(B,1,ecg_h);
% fvtool(B)

%% Create and plot power spectra of entire ecg_h
L = length(ecg_h);
NFFT = 2^nextpow2(L);
Y = fft(ecg_h,NFFT) / L;
Y = abs(fftshift(Y));
f = F_s * linspace(-0.5,0.5,NFFT);

figure;
plot(f,Y)
axis([0 100 0 1.2 * max(Y)]);
hold on;

title('Amplitude Spectrum of ecg_h')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

xval = (1 / meanDuration)/2;
x=[xval,xval];
y=[0,1000];
plot(x,y,'-r')
plot(-x,y,'-r')

%% Visualize series of heartbeats vertically
heartbeats = 3;
offset = 0;
figure;
for i=1:(heartbeats)
    subplot(heartbeats,1,i);
    nrange = (t_h(i + offset):t_h(i + offset + 1));
    plot(T*nrange,ecg_h(nrange));
    title(num2str(i+offset));
    axis([T*nrange(1) T*nrange(end) -1000 1000]);
end

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
    maxlen = max(length(heartbeat1),length(heartbeat2));
    
    heartbeat1 = interp1(heartbeat1,linspace(1,numel(heartbeat1),maxlen));
    heartbeat2 = interp1(heartbeat2,linspace(1,numel(heartbeat2),maxlen));
    
    subplot(3,pairs,i);
    plot(T * (1:length(heartbeat1)),heartbeat1);
    axis([0 maxlen*T -inf inf]);
    title(num2str(k+offset));
    xlabel('time /s')
    ylabel('amplitude /mV')
    hold on;
    xval = 0.174;
    x=[xval,xval];
    y=[-100,1000];
    plot(x,y,'-r')
    plot(-x,y,'-r')
    
    subplot(3,pairs,i + pairs);
    plot(T * (1:length(heartbeat2)),heartbeat2);
    axis([0 maxlen*T -inf inf]);
    title(num2str(k+1+offset));
    xlabel('time /s')
    ylabel('amplitude /mV')
    hold on;
    xval = 0.174
    x=[xval,xval];
    y=[-100,1000];
    plot(x,y,'-r')
    plot(-x,y,'-r')
    
    diff = abs(heartbeat1 - heartbeat2);
    subplot(3,pairs,i + 2 * pairs);
    plot(T * (1:length(diff)),diff);
    axis([0 maxlen*T -inf inf]);
    title(['|' num2str(k+offset) ' - ' num2str(k+offset + 1) '|']);
    xlabel('time /s')
    ylabel('amplitude /mV')
    hold on;
    xval = 0.174
    x=[xval,xval];
    y=[-100,1000];
    plot(x,y,'-r')
    plot(-x,y,'-r')
end
