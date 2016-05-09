F_s = 1000;
twave_length = 300; %samples
twave_freq = 2.5; %Hz
std_twave = sin(2 * pi * (twave_freq/F_s) * (0:(twave_length - 1))) + 1;

N = 500; %heartbeats
STT_matrix = zeros(N,twave_length);
TWA_episodes = [50:100 200:210 300:400];
TWA_amplitude = 1;
for n = 1:N
    if any(n==TWA_episodes)
        STT_matrix(n,:) = std_twave + (TWA_amplitude * (-1)^n);
    else
        STT_matrix(n,:) = std_twave;
    end
end
surf(STT_matrix)