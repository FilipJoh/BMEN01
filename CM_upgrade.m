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
threshold=1;
%% Pig one
% ecg_p=ecg_p1;
% t=t_p1;
% nbrOfBeats=length(t)-2;
% threshold=1e-3;

%% Pig two
% ecg_p=ecg_p2;
% t=t_p2;
% nbrOfBeats=length(t)-2; 
% threshold=0.01;

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
tempSignalMat=zeros(128,TwindowLength);
signalMat=zeros(nbrOfBeats,TwindowLength);%zeros(nbrOfBeats,twave_end-twave_beg);
heartBeatMatrix=cell(nbrOfBeats,1);
starts=zeros(128,1);
%figure;
%0.23
figure;
for i=1:size(ecg_p,1)
    ecg=ecg_p(i,:);
    for a=1:nbrOfBeats
        heartBeatMatrix{a}=ecg(t(a):t(a+1)); 
    end
    
    for a = 1:(nbrOfBeats-127)
        for v=1:128
            if beatDuration(a)*T<0.6
                startoffset=0.06;
            elseif beatDuration(a)*T<1.1
                startoffset=0.1;    
            else
                startoffset=0.15;
            end
                    
            twave_beg =(startoffset / T);
            %twave_end = twave_beg + ((0.35-startoffset) / T);
            %temp = ecg(twave_beg:twave_end);
            tempSignalMat(v,:)=heartBeatMatrix{a+v-1}(twave_beg:(twave_beg-1+TwindowLength));
            starts(v)=twave_beg;
        end
        TempMed=median(tempSignalMat);
        for v=1:128
            Best=0;
            bestStart=starts(v);
            for u=-0.03:0.001:0.03
                offsetVal=TempMed*(heartBeatMatrix{v}((starts(v)+floor(u/T)):(starts(v)+floor(u/T)-1+TwindowLength)))';
                if offsetVal>Best;
                    Best=offsetVal;
                    bestStart=starts(v)+u/T;
                end    
            end
            signalMat(v+a-1,:)=heartBeatMatrix{v+a-1}(bestStart:(bestStart-1+TwindowLength));
        end
    end
% create bandstopfilter in order to sort of respiratory signal
[n,Wn]=buttord([0.11 0.38],[0.14 0.35],3,60);
[B,A]=butter(n,Wn,'stop');
% fvtool(B,A);

for a=1:size(nbrOfBeats)
    y=filtfilt(B,A,signalMat(a,:)); 
end
    
    
%detextion    
    window=7;
%     ACI=zeros(window,1);
%     ACM=zeros(window,1);
    
    ACI=zeros(nbrOfBeats,1);
    ACM=zeros(nbrOfBeats,1);
    ACMmatrix=[];
    startWin=1;
    endWin=startWin+window-1;
    alt=zeros(size(signalMat,1),1);
    %figure;
    while endWin~=size(signalMat,1);
        Tmdn=median(signalMat(startWin:endWin,:));
        for j=startWin:endWin%1:window     
               % ACI(j)=max(xcorr(signalMat(startWin-1+j,:),Tmdn))./max(xcorr(Tmdn,Tmdn));%signalMat(a,:)*Tmdn'./(Tmdn*Tmdn');
                %ACM(j)=2*abs(ACI(j)-1)*sum(Tmdn.^2)./sum(abs(Tmdn));
                [acf,lags,bound]=autocorr(Tmdn);
                ACI(j)=max(xcorr(signalMat(j,:),Tmdn))/max(xcorr(Tmdn,Tmdn));
                %ACI(j)=sum(signalMat(j,:).*Tmdn)/(sum(Tmdn.^2));%max(xcorr(signalMat(j,:),Tmdn))./max(xcorr(Tmdn,Tmdn));%signalMat(a,:)*Tmdn'./(Tmdn*Tmdn');
                ACM(j)=2*abs(ACI(j)-1)*sum(Tmdn.^2)./sum(abs(Tmdn));
        end
        ACMmatrix=[ACMmatrix ACM];
%         plot(1:length(ACI),ACI-1)
%         axis([0 nbrOfBeats -1 1])
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
    plot(1:length(alt),alt,'LineWidth',2);
    for j=1:size(ACM,1)
        if alt(j)==1
            ACMcomp(j)=mean(ACMmatrix(j,find(ACMmatrix(j,:))));
        end
    end
    hold on;
    plot(1:length(ACMcomp),ACMcomp,'LineWidth',1)
    %figure;
    %contourf(ACMmatrix,'EdgeColor','none')
    title(['ECG-lead no: ' num2str(i)])
    xlabel('heartbeats')
    ylabel('TWA amplitude (mv)')
%     figure;
%     plot((1:length(temp))*T,temp);

  
    % figure;
    % plot(1:length(Tmdn),Tmdn)



    %[B,A]=butter(16,0.5,'high');
    %ACI=filter(B,A,ACI);

    
    
%     ACM=medfilt1(ACM,30);

%     figure;
%     if(size(ecg_p,1)>1)
%            subplot(size(ecg_p,1)/2,2,i)
%     end
end