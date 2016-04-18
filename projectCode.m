%% load data
load('AlternansData.mat');
T=1/1000*(1:length(ecg_h));
heartbeats=6;
heartbeats=heartbeats+1;
range=t_h(1):(t_h(heartbeats)-1);

freqD=fft(ecg_h);
figure;
freqD=abs(ecg_h);
plot(freqD);

figure;
plot(T(range),ecg_h(range))
axis([T(range(1)) T(range(end)) -1000 1000])

figure;
hold on;
legendstring={};
col=hsv(heartbeats-1);
for i=2:heartbeats
    %figure;
    rangeb=t_h(i-1):(t_h(i)-1);
    ecgRange=ecg_h(rangeb);
    normalizedRange=rangeb./rangeb(end);
    plot(normalizedRange,ecgRange,'color',col(i-1,:))
    axis([0 1 0 700])
    legendstring{i-1}=strcat('heartbeat no: ',num2str(i-1));
end
legend(legendstring)