clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
signal = ecg_h;
t_beat = t_h;

df = 2;
signal = decimate(signal,df);
t_beat = round(t_beat / df);
F_s = F_s / df;

cutFreq = 50; %Hz
order = 100;
B = fir1(order,cutFreq/(F_s/2));
signal = filtfilt(B,1,signal);

% Singlebin
twave_beg = 0.06; %s
twave_end = 0.29; %s
[alt_amp,~] = CD_singlebin(signal,t_beat,twave_beg,twave_end,F_s);
figure;
plot(1:length(alt_amp),alt_amp)
title('Alternans amplitude');
xlabel('Heartbeat')

% Multibin
twave_beg = 0.06; %s
twave_end = 0.4; %s
binDuration = 0.01; %s
[alt_amp,~,tih,n] = CD_multibin(signal,t_beat,binDuration,twave_beg,twave_end,F_s);
n = (n./60);
figure;
contourf(tih,n,alt_amp)
xlabel('Time after R wave / ms')
ylabel('Time / minutes')
c=colorbar;
xlabel(c,'TWA amplitude / \muV');