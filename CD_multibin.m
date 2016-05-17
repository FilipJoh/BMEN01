function [alt_amp,alt_phase,t,n] = CD_multibin(signal,t_beat,binDuration,twave_beg,twave_end,F_s)
%CD Summary of this function goes here
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
    
    mdnMat(n,:) = twave;
end
twave_mdn = median(mdnMat);

% Generate main time series
binSize = binDuration / T; %samples
nbrOfBins = floor(((twave_end - twave_beg) / binDuration) + 1e-10);
X = zeros(N,nbrOfBins);
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
    
    for k = 1:nbrOfBins
        bin_samples = (1 + (k - 1)*binSize):((k)*binSize);
        bin_amplitudes = twave(bin_samples);
        bin_times = twave_beg + T.*(bin_samples-1);
        I = cumtrapz(bin_times,bin_amplitudes);
        X(n,k) = I(end) / (binDuration);
    end
end

%remove low frequency components
[B,A] = butter(16,0.9,'high');
X = (filtfilt(B,A,X));

% Complex demodulation
altfreq = 0.5;
n = 0:N-1;
complex_sin = exp(1i*2*pi*altfreq * n)';
complex_sin_mat = repmat(complex_sin,1,nbrOfBins);
alternans = X .* complex_sin_mat;

% keep only low frequency component
[B,A] = butter(16,0.1,'low');
alternans = filtfilt(B,A,alternans);

alt_amp = abs(alternans);
alt_phase = angle(alternans);
t = twave_beg + binDuration*(0:(nbrOfBins - 1));

n = T .* t_beat(1:N);
end
