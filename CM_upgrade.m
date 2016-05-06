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

% %% Human
% ecg_p=ecg_h;
% t=t_h;
% nbrOfBeats=find(t<length(ecg_p),1,'last')-1;
% threshold=1;
% %% Pig one
% ecg_p=ecg_p1;
% t=t_p1;
% nbrOfBeats=length(t)-2;
% threshold=1e-3;

%% Pig two
ecg_p=ecg_p2;
t=t_p2;
nbrOfBeats=length(t)-2; 
threshold=0.01;

% create a lowpass filter
cutFreq = 50; %Hz
B = fir1(100,cutFreq/(F_s/2));
[H,w]=freqz(B,1);
%plot(F_s * w/2/pi,abs(H))
ecg_p = filtfilt(B,1,ecg_p')';

%% calculations

 %698;%length(t_h)-3;
beatDuration=t(2:nbrOfBeats+1)-t(1:nbrOfBeats);
maxBeatDuration=max(beatDuration)+1;
twave_beg = t(1) + (0.04 / T);
twave_end = twave_beg + (0.35 / T);

signalMat=zeros(nbrOfBeats,twave_end-twave_beg);
%figure;
%0.23
figure;
for i=1:size(ecg_p,1)
    ecg=ecg_p(i,:);
    for a = 1:nbrOfBeats
            twave_beg = t(a) + (0.06 / T);
            twave_end = twave_beg + ((0.35-0.06) / T);
            temp = ecg(twave_beg:twave_end);
            %temp=ecg_h(t_h(a):t_h(a+1));
            %plot(1:length(temp),temp);
            %temp=temp(0.06/T:0.23/T);
            %plot(1:length(temp),temp);

            %if length(temp) < maxBeatDuration
%                 signalMat(a,1:maxBeatDuration) = interp1(temp,linspace(1,numel(temp),maxBeatDuration));
%             else   
                signalMat(a,1:length(temp))=temp;
%             end
%             plot((1:length(temp))*T,temp);
    end
    
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
    ylabel('Alterans (On/Off)')
%     figure;
%     plot((1:length(temp))*T,temp);

  
    % figure;
    % plot(1:length(Tmdn),Tmdn)



    [B,A]=butter(16,0.5,'high');
    ACI=filter(B,A,ACI);

    
    
%     ACM=medfilt1(ACM,30);

%     figure;
%     if(size(ecg_p,1)>1)
%            subplot(size(ecg_p,1)/2,2,i)
%     end
end