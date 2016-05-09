function [t,TWA] = simulateTwave(T,AlternansRangeTime )
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
Alternans=ones(1,length(samples));
AlternansFreq=TwaveFreq/2;
AlternansAmp=1;
AlternansOffset=1;
AlternansRangeSamples=floor(AlternansRangeTime./T);
Alternans(AlternansRangeSamples(1):AlternansRangeSamples(2))=sin(2*pi*AlternansFreq*samples(AlternansRangeSamples(1):AlternansRangeSamples(2)))+AlternansOffset;
TWA=Twave.*Alternans;
t=1:beatSamples:(128*beatSamples);
end

