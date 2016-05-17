function [alt_amp,alt_phase] = CD_singlebin(signal,t_beat,twave_beg,twave_end,F_s)
%CD_SINGLEBIN Summary of this function goes here
%   Detailed explanation goes here
T = 1 / F_s;
    
% Generate median t-wave
N = find(t_beat < length(signal),1,'last') - 1; %length(t_beat) exceeds signal length?
mdnMat = zeros(N,round((twave_end - twave_beg) * F_s));
for n = 1:N
    heartbeat = signal(t_beat(n):(t_beat(n + 1) - 1));
    twave_beg_index = 1 + (twave_beg / T);
    twave_end_index = twave_beg_index + ((twave_end - twave_beg) / T) - 1;
    twave = heartbeat(twave_beg_index:twave_end_index);
    
    baseline = 0; %todo
    twave = twave - baseline;
    
%     cutFreq = 20; %Hz
%     order = 16;
%     [B,A] = butter(order,cutFreq/(F_s/2));
%     twave = filtfilt(B,A,twave);
    
    mdnMat(n,:) = twave;
end
twave_mdn = median(mdnMat);

% Generate main time series
X = zeros(N,1);
for n = 1:N
    heartbeat = signal(t_beat(n):(t_beat(n + 1) - 1));
    twave_beg_index = 1 + (twave_beg / T);
    twave_end_index = twave_beg_index + ((twave_end - twave_beg) / T) - 1;
    
    % t-wave alignment
    corr_max = 0;
    best_lag = 0;
    for lag = (-0.03*F_s):(0.03*F_s)
        corr_curr = twave_mdn * heartbeat((twave_beg_index + lag):(twave_end_index + lag))';
        if corr_curr > corr_max
            corr_max = corr_curr;
            best_lag = lag;
        end
    end
    twave = heartbeat((twave_beg_index + best_lag):(twave_end_index + best_lag));
    
    baseline = 0; %todo
    twave = twave - baseline;
    
    cutFreq = 20; %Hz
    order = 16;
    [B,A] = butter(order,cutFreq/(F_s/2));
    twave = filtfilt(B,A,twave);
    
    I = cumtrapz(twave_beg:T:(twave_end-T),twave);
    X(n) = I(end) / (T * length(twave));   
end

%remove low frequency components
[B,A] = butter(16,0.9,'high');
X = filtfilt(B,A,X);

% Complex demodulation
altfreq = 0.5;
n = 0:N-1;
alternans = X .* 2 .* exp(1i*2*pi*altfreq * n)';

% keep only low frequency component
[B,A] = butter(16,0.1,'low');
alternans = filtfilt(B,A,alternans);

alt_amp = abs(alternans);
alt_phase = angle(alternans);
end

