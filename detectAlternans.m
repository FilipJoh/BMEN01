function [ boolVec ] = detectAlternans( ACI,threshold )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
counter=0;
% threshold=0.2;
boolVec=zeros(length(ACI),1);
detect=7;
for i=3:length(ACI)
     if(ACI(i-1)-ACI(i-2)>threshold && ACI(i)-ACI(i-1)< -threshold) || (ACI(i-1)-ACI(i-2)<-threshold && ACI(i)-ACI(i-1)>threshold)
         counter=counter+1;
     else
         counter=0;
     end
     if counter==detect
        boolVec(i-(detect+1):i)=1; 
     elseif counter >detect
        boolVec(i)=1;
     end
end

end

