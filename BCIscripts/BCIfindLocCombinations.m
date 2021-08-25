function sAsVconds = BCIfindLocCombinations(trueLocsAV)
% Find all unique combinations of sA and sV

%Separate auditory and visual locations
sA_all = trueLocsAV(:,1);
sV_all = trueLocsAV(:,2);

%Find number of unique bisensory conditions (combinations of sA and sV)   
sA_unique = unique(sA_all(~isnan(sA_all)));
sV_unique = unique(sV_all(~isnan(sV_all)));
sAsVconds = [];
for i=1:numel(sA_unique)
    i_A = sA_all == sA_unique(i);
    for j=1:numel(sV_unique)
        i_V = sV_all == sV_unique(j);
        if sum(i_A & i_V) > 0
            sAsVconds = [sAsVconds; sA_unique(i) sV_unique(j)];
        end
    end
end

%Add unisensory visual conditions (A == NaN)
if sum(isnan(sA_all)) > 0
    i_A = isnan(sA_all);
    for j=1:numel(sV_unique)
        i_V = sV_all == sV_unique(j);
        if sum(i_A & i_V) > 0
            sAsVconds = [sAsVconds; NaN sV_unique(j)];
        end
    end
end

%Add unisensory auditory conditions (V == NaN)
if sum(isnan(sV_all)) > 0
    i_V = isnan(sV_all);
    for j=1:numel(sA_unique)
        i_A = sA_all == sA_unique(j);
        if sum(i_A & i_V) > 0
            sAsVconds = [sAsVconds; sA_unique(j) NaN];
        end
    end
end

end %[EOF]