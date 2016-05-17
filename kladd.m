clc;
clear all;
% close all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;
leads = ecg_p1;
t_beat = t_p1;

%% Ahuba
signal = leads(9,:);
pairs = 5;
offset = 315;
figure;
for i=1:(pairs)
    k = 2*i - 1;
    range1 = (t_beat(k + offset):t_beat(k + offset + 1));
    range2 = (t_beat(k + 1 + offset):t_beat(k + 1 + offset + 1));
    heartbeat1 = signal(range1);
    heartbeat2 = signal(range2);
    maxlen = max(length(heartbeat1),length(heartbeat2));
    
    heartbeat1 = interp1(heartbeat1,linspace(1,numel(heartbeat1),maxlen));
    heartbeat2 = interp1(heartbeat2,linspace(1,numel(heartbeat2),maxlen));
    
    subplot(3,pairs,i);
    plot(T * (1:length(heartbeat1)),heartbeat1);
    axis([0 maxlen*T -inf inf]);
    title(num2str(k+offset));
    xlabel('time /s')
    ylabel('amplitude /mV')
    
    subplot(3,pairs,i + pairs);
    plot(T * (1:length(heartbeat2)),heartbeat2);
    axis([0 maxlen*T -inf inf]);
    title(num2str(k+1+offset));
    xlabel('time /s')
    ylabel('amplitude /mV')
    
    diff = abs(heartbeat1 - heartbeat2);
    subplot(3,pairs,i + 2 * pairs);
    plot(T * (1:length(diff)),diff);
    axis([0 maxlen*T -inf inf]);
    title(['|' num2str(k+offset) ' - ' num2str(k+offset + 1) '|']);
    xlabel('time /s')
    ylabel('amplitude /mV')
end


%% Heartrate
T = 1 / F_s;
durationsSamples = (t_beat(2:end) - t_beat(1:(end-1)));
durations = T * durationsSamples;

figure;
plot(1:length(durations),durations);
xlabel('heartbeat')
ylabel('duration /s')

meanDuration = mean(durations);
meanHeartrate = 60*(1 / meanDuration)

