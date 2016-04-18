%% load data
load('AlternansData.mat');
T=1/1000*(1:length(ecg_h));
heartbeats=6;
heartbeats=heartbeats+1;
range=t_h(1):(t_h(heartbeats)-1);

figure;
plot(T(range),ecg_h(range))
axis([T(range(1)) T(range(end)) -1000 1000])

figure;
hold on;
for i=2:heartbeats
    rangeb=[t_h(heartbeats-1) t_h(heartbeats)-1];
    plot()
end