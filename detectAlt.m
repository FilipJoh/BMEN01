function [ alternans ] = detectAlt(ACI)
    if min(ACI(1:2:end) >= 0) & min(ACI(2:2:end-1)< 0) || min(ACI(1:2:end)) >= 0 & min(ACI(2:2:end-1)) < 0
        alternans=ones(length(ACI),1);
    else    
        alternans=zeros(length(ACI),1);
    end
end

