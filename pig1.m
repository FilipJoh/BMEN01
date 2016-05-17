%% Initialize
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
leads = ecg_p1;
t = t_p1;
t_beat = t;

%% Find alternans with CD

% Singlebin
figure;
for i = 1:12
   signal = leads(i,:);
   t_beat = t;
   F_s = 1000;
   twave_beg = 0.1; %s
   twave_end = 0.35; %s
   
   df = 2;
   signal = decimate(signal,df);
   t_beat = round(t_beat / df);
   F_s = F_s / df;
   
   cutFreq = 50; %Hz
   order = 100;
   B = fir1(order,cutFreq/(F_s/2));
   signal = filtfilt(B,1,signal);
   
   [alt_amp,alt_phase] = CD_singlebin(signal,t_beat,twave_beg,twave_end,F_s);
   subplot(3,4,i);
   plot(1:length(alt_amp),alt_amp)
   title(['Alternans lead ' num2str(i)]);
   xlabel('Heartbeat')
end

% Multibin
figure;
for i = 1:12
   signal = leads(i,:);
   t_beat = t;
   F_s = 1000;
   twave_beg = 0.06; %s
   twave_end = 0.4; %s
   binDuration = 0.01; %s
   
   df = 2;
   signal = decimate(signal,df);
   t_beat = round(t_beat / df);
   F_s = F_s / df;
   
   cutFreq = 50; %Hz
   order = 100;
   B = fir1(order,cutFreq/(F_s/2));
   signal = filtfilt(B,1,signal);
   
   [alt_amp,alt_phase,tih,n] = CD_multibin(signal,t_beat,binDuration,twave_beg,twave_end,F_s);
   n = (n./60);
   subplot(3,4,i);
   contourf(tih,n,alt_amp)
   title(['ECG lead no: ' num2str(i)]);
   xlabel('Time after R wave / ms')
   ylabel('Time / minutes')
   c=colorbar;
   xlabel(c,'TWA amplitude / \muV');
end
