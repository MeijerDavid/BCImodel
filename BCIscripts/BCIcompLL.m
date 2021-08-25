function LLvector = BCIcompLL(paramsConv,P,trueLocsAV,responsesAVC,i_Conditions,refLL)
% This function returns a row vector with log likelihoods for each response

%Check refLL input parameter (serves as a break from further computing LLs - used during MCMC)   
if nargin < 6
    refLL = [];
end
refLL_break = false;   %Initialize the boolean

%Convert the parameters back to 'normal' space
paramsReal = BCIconvertParams2Fit(paramsConv,P,'conv2real');

%Gather all unique locations (combinations of sA and sV)
LocConds = P.sAsVconds;
nLocConds = size(LocConds,1);

%Initialize a matrix wherein we save the LL of each response
LLresponses = nan(size(responsesAVC));

%Loop across all experimental conditions
for c=1:P.nConditions
    
    %Overwrite the parameters to fit for this condition
    Pcopy = P;
    for i=1:numel(P.ParamsPerCond{1,c})
        paramIdx = P.ParamsPerCond{1,c}(i);
        Pcopy.(P.ParamNames2Fit{1,paramIdx}) = paramsReal(paramIdx);
    end

    %Loop across all locations (combinations of sA and sV)
    for i=1:nLocConds

        %Select the true locations of the stimuli
        sA = LocConds(i,1);
        sV = LocConds(i,2);
        
        %Select the relevant observer's responses
        if isnan(sA) 
            i_AVcond = i_Conditions(:,c) & isnan(trueLocsAV(:,1)) & (trueLocsAV(:,2) == sV);    %unisensory visual
        elseif isnan(sV)
            i_AVcond = i_Conditions(:,c) & (trueLocsAV(:,1) == sA) & isnan(trueLocsAV(:,2));    %unisensory auditory
        else
            i_AVcond = i_Conditions(:,c) & (trueLocsAV(:,1) == sA) & (trueLocsAV(:,2) == sV);   %audiovisual
        end
        
        %If there are trials in this condition+location combination
        if sum(i_AVcond) > 0
            
            %Select the relevant observer's responses for this condition
            i_relResp_A = i_AVcond & ~isnan(responsesAVC(:,1));
            i_relResp_V = i_AVcond & ~isnan(responsesAVC(:,2));
            i_relResp_C = i_AVcond & ~isnan(responsesAVC(:,3));
            sA_resp = responsesAVC(i_relResp_A,1);                          %These can be 'empty'
            sV_resp = responsesAVC(i_relResp_V,2);
            CSJ_resp = responsesAVC(i_relResp_C,3);
            
            %Compute log-likelihoods of the responses
            if sum([numel(sA_resp); numel(sV_resp); numel(CSJ_resp)]) > 0   %If there are valid responses...
                [LLresponses(i_relResp_A,1), LLresponses(i_relResp_V,2), LLresponses(i_relResp_C,3)] = BCIcompLL_bci(Pcopy,sA,sV,sA_resp,sV_resp,CSJ_resp);
            end
        end
        
        %Check if already smaller than reference LL - if so, no need to further compute LLs (small speed-up for use during MCMC)
        if ~isempty(refLL)
            refLL_break = (nansum(LLresponses,'all') < refLL);
            if refLL_break
                break; %break from locations for-loop
            end
        end
    end
    
    %Also break out of conditions loop if already smaller than reference LL 
    if refLL_break 
        break; 
    end
    
end

%Create one long row vector with log-likelihoods for all responses (i.e. ignoring NaNs)   
LLvector = [LLresponses(~isnan(LLresponses(:,1)),1); LLresponses(~isnan(LLresponses(:,2)),2); LLresponses(~isnan(LLresponses(:,3)),3)]';

%Avoid returning -inf as a probability: e.g. coming from LL = log(prob = 0) = -inf. 
%This is not a good way of avoiding it. Instead, it would be better to adjust the parameter bounds.
%LLvector(LLvector == -inf) = -realmax;

end %[EOF]
