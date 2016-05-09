
clc;
clear all;
F_s = 1000;
T = 1 / F_s;

%create ordinary twave
%Twave duration
    TwaveDur=200e-3;
%Twave frequency
    TwaveFreq=1/TwaveDur;
%Samples per beat 
beatSamples=TwaveDur/T;

%total samples
samples=1:T:(132*TwaveDur);
Twave=-sin(2*pi*TwaveFreq*samples+1);


%create alternans
A=ones(1,length(samples));
AlternansFreq=1/2;
AlternansAmp=0.5;
AlternansOffset=1;
%AlternansRangeSamples=floor([5.0 10.0]./T);
%A(AlternansRangeSamples(1):AlternansRangeSamples(2))=AlternansAmp*(-1).^(AlternansRangeSamples(1):AlternansRangeSamples(2))+AlternansOffset;
TWA=Twave.*A;
t=1:beatSamples:(128*beatSamples);
nbrOfBeats=length(t)-2;


signalMat=zeros(127,beatSamples+1);
for a=1:127
   signalMat(a,:)=TWA(t(1):t(2));  
   if a>40 && a<80
        signalMat(a,:)=(1*(-1)^a)+signalMat(a,:);
   end    
end
ecg_p=TWA;

% Tmdn=median(signalMat);
 window=7;   
    ACI=zeros(nbrOfBeats,1);
    ACM=zeros(nbrOfBeats,1);
    ACMmatrix=[];
    ACImatrix=[];
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
        ACImatrix=[ACImatrix ACI];
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
    ACIsamples=find(alt);
    ACItime=t(find(alt));
    time=t(1:length(ACMcomp));
    
    for j=1:size(ACM,1)
        if alt(j)==1
            ACMcomp(j)=mean(ACMmatrix(j,find(ACImatrix(j,:))));
        end
    end
    hold on;
    plot(samples,ACMcomp,'r','LineWidth',1)
    hold on;
    plot(ACIsamples,zeros(length(ACItime),1),'b.','LineWidth',5);
    
    %figure;
    %contourf(ACMmatrix,'EdgeColor','none')
    title(['ECG lead no: ' num2str(i)])
    xlabel('Time /minutes')
    ylabel('TWA amplitude  /\muv')
    legend('TWA estimate','Binary mask')