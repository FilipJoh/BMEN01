%% load data and define constants
clc;
clear all;
load AlternansData.mat; % ecg_h , ecg_p1, ecg_p2, t_h, t_p1, t_p2
F_s = 1000;
T = 1 / F_s;

%% create a lowpass filter
cutFreq = 50; %Hz
B = fir1(100,cutFreq/(F_s/2));
[H,w]=freqz(B,1);
%plot(F_s * w/2/pi,abs(H))
ecg_h = filtfilt(B,1,ecg_h);



%% TEST

%% Human
ecg_p=ecg_h;
t=t_h;
nbrOfBeats=find(t<length(ecg_p),1,'last')-1;
threshold=0.1;
%% Pig one
% ecg_p=ecg_p1;
% t=t_p1;
% nbrOfBeats=length(t)-2;
% threshold=1e-3;

%% Pig two
% ecg_p=ecg_p2(1,:);
% t=t_p2;
% nbrOfBeats=length(t)-2; 
% threshold=0.01;

% create a lowpass filter
% cutFreq = 50; %Hz
% B = fir1(100,cutFreq/(F_s/2));
% [H,w]=freqz(B,1);
% plot(F_s * w/2/pi,abs(H))
% ecg_p = filtfilt(B,1,ecg_p')';

%% calculations

 %698;%length(t_h)-3;
beatDuration=t(2:nbrOfBeats+1)-t(1:nbrOfBeats);
maxBeatDuration=max(beatDuration)+1;
twave_beg = t(1) + (0.04 / T);
twave_end = twave_beg + (0.35 / T);

signalMat=zeros(nbrOfBeats,twave_end-twave_beg);
figure;
%0.23

for i=1:size(ecg_p,1)
    ecg=ecg_p(i,:);
    for a = 1:nbrOfBeats
            twave_beg = t(a) + (0.06 / T);
            twave_end = twave_beg + (0.35 / T);
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
%     figure;
%     plot((1:length(temp))*T,temp);

    Tmdn=median(signalMat);
    % figure;
    % plot(1:length(Tmdn),Tmdn)


    ACI=zeros(nbrOfBeats,1);
    for a=1:nbrOfBeats
         ACI(a)=max(xcorr(signalMat(a,:),Tmdn))./max(xcorr(Tmdn,Tmdn));%signalMat(a,:)*Tmdn'./(Tmdn*Tmdn');
         ACM(a)=2*abs(ACI(a)-1)*sum(Tmdn.^2)./sum(abs(Tmdn));
    end

    [B,A]=butter(16,0.5,'high');
    ACI=filter(B,A,ACI);

    boolVec=detectAlternans(ACI,threshold);
    
%     ACM=medfilt1(ACM,30);

%     figure;
    if(size(ecg_p,1)>1)
%           subplot(size(ecg_p,1)/2,2,i)
    end
 
    plot(1:length(ACM),ACM);
    %hold on;
    %figure;
    %plot(1:length(Tmdn),signalMat(498,:),'b',1:length(Tmdn),Tmdn,'r');
    %figure;
    %plot(1:length(ACM),ACM);
%     plot(1:length(ACM),boolVec)
    %legend(['AutoCorrelation Index'])
    title(['ACM vs. time lead: ' num2str(i)]);
   % hold off;
end



