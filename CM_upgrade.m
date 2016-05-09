%% load data and define constants
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;

%% create a lowpass filter
% cutFreq = 50; %Hz
% B = fir1(100,cutFreq/(F_s/2));
% [H,w]=freqz(B,1);
% %plot(F_s * w/2/pi,abs(H))
% ecg_h = filtfilt(B,1,ecg_h);
% 


%% TEST

%% Human
ecg_p=ecg_h;
t=t_h;
nbrOfBeats=find(t<length(ecg_p),1,'last')-1;

%% Pig one
% ecg_p=ecg_p1;
% t=t_p1;
% nbrOfBeats=length(t)-2;


%% Pig two
% ecg_p=ecg_p2;
% t=t_p2;
% nbrOfBeats=length(t)-2; 

% create a lowpass filter
cutFreq = 50; %Hz
B = fir1(100,cutFreq/(F_s/2));
[H,w]=freqz(B,1);
%plot(F_s * w/2/pi,abs(H))
ecg_p = filtfilt(B,1,ecg_p')';

%% calculations
% Preprocessing

 %698;%length(t_h)-3;
beatDuration=t(2:nbrOfBeats+1)-t(1:nbrOfBeats);
meanDuration=mean(beatDuration)*T;
maxBeatDuration=max(beatDuration)+1;
twave_beg = t(1) + (0.04 / T);
twave_end = twave_beg + (0.35 / T);
TwindowLength=floor((0.4*sqrt(meanDuration))/T);
tempSignalMat=zeros(nbrOfBeats,TwindowLength);
signalMat=zeros(nbrOfBeats,TwindowLength);%zeros(nbrOfBeats,twave_end-twave_beg);
heartBeatCell=cell(nbrOfBeats,1);
starts=zeros(nbrOfBeats,1);

figure;
for i=1:size(ecg_p,1)
    ecg=ecg_p(i,:);
    for a=1:nbrOfBeats
        heartBeatCell{a}=ecg(t(a):t(a+1)); 
    end
    
    for v=1:nbrOfBeats
        if beatDuration(a)*T<0.6
            startoffset=0.06;
        elseif beatDuration(a)*T<1.1
            startoffset=0.1;    
        else
            startoffset=0.15;
        end

        twave_beg =(startoffset / T);
        tempSignalMat(v,:)=heartBeatCell{v}(twave_beg:(twave_beg-1+TwindowLength));
        starts(v)=twave_beg;
    end
    TempMed=median(tempSignalMat);
    for v=1:nbrOfBeats
        Best=0;
        bestStart=starts(v);
        for u=-0.03/T:0.03/T
            offsetVal=TempMed*(heartBeatCell{v}((starts(v)+floor(u)):(starts(v)+floor(u)-1+TwindowLength)))';
            if offsetVal>Best;
                Best=offsetVal;
                bestStart=starts(v)+u;
            end    
        end
        signalMat(v,:)=heartBeatCell{v}(bestStart:(bestStart-1+TwindowLength));
    end

[n,Wn]=buttord([0.1 0.39],[0.14 0.35],3,60);
[B,A]=butter(16,Wn,'stop');

for a=1:size(nbrOfBeats)
    %y=filtfilt(B,A,signalMat(a,:)); 
end
    
   
%detection    
    window=15;   
    ACI=zeros(nbrOfBeats,1);
    ACM=zeros(nbrOfBeats,1);
    ACMmatrix=[];
    startWin=1;
    endWin=startWin+window-1;
    alt=zeros(size(signalMat,1),1);

    while endWin~=size(signalMat,1);
        Tmdn=median(signalMat(startWin:endWin,:));
        for j=startWin:endWin     
                [acf,lags,bound]=autocorr(Tmdn);
                ACI(j)=max(xcorr(signalMat(j,:),Tmdn))/max(xcorr(Tmdn,Tmdn)); %ACI(j)=sum(signalMat(j,:).*Tmdn)/(sum(Tmdn.^2));%max(xcorr(signalMat(j,:),Tmdn))./max(xcorr(Tmdn,Tmdn));%signalMat(a,:)*Tmdn'./(Tmdn*Tmdn');
                ACM(j)=2*abs(ACI(j)-1)*sum(Tmdn.^2)./sum(abs(Tmdn));
        end
        ACMmatrix=[ACMmatrix ACM];
        alt(startWin:endWin)=or(alt(startWin:endWin),detectAlt(ACI(startWin:endWin)-1));%mean(ACM(startWin:endWin));
        ACM=zeros(nbrOfBeats,1);
        ACI=zeros(nbrOfBeats,1);
        startWin=startWin+1;
        endWin=endWin+1;
        
        
    end
    %figure;
    if(size(ecg_p,1)>1)
           subplot(3,4,i)
    end
    ACMcomp=zeros(length(ACI),1);
    
 
    %figure;
    samples=1:length(ACMcomp);
    ACItime=t(find(alt))*T/60;
    time=t(1:length(ACMcomp))*T/60; %minutes
    
    for j=1:size(ACM,1)
        if alt(j)==1
            ACMcomp(j)=mean(ACMmatrix(j,find(ACMmatrix(j,:))));
        end
    end
    hold on;
    plot(time,ACMcomp,'r','LineWidth',1)
    hold on;
    plot(ACItime,zeros(length(ACItime),1),'b.','LineWidth',5);
    
    %figure;
    %contourf(ACMmatrix,'EdgeColor','none')
    title(['ECG lead no: ' num2str(i)])
    xlabel('Time /minutes')
    ylabel('TWA amplitude  /\muv')
    legend('TWA estimate','Binary mask')
end