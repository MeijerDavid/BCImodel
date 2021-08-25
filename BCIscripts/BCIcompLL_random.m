function [LLvectorA,LLvectorV,LLvectorC] = BCIcompLL_random(P,sA,sV,sA_resp,sV_resp,CSJ_resp)
% Generate random response frequencies for one condition (combination of sA and sV)
% And compute loglike for A, V and CSJ responses

% Gather relevant variables from structure P
RespLocs = P.RespLocs;
RespRange = P.RespRange;

% Determine number of trials
nTrialsA = numel(sA_resp);
nTrialsV = numel(sV_resp);
nTrialsC = numel(CSJ_resp);

% Initialize output vectors
LLvectorA = nan(nTrialsA,1);
LLvectorV = nan(nTrialsV,1);
LLvectorC = nan(nTrialsC,1);

%Discretized responses
if all(~isnan(RespLocs))
    
    %All response buttons have equal probability
    nRespLocs = numel(RespLocs);
    likeFunAV = (1/nRespLocs)*ones(1,nRespLocs);
    if (nTrialsA > 0); LLvectorA(:) = log(likeFunAV(sA_resp)); end          %Note that the responses were already transformed to bin_idx numbers before
    if (nTrialsV > 0); LLvectorV(:) = log(likeFunAV(sV_resp)); end
    
%Continuous responses  
else
    
    %Compute the constant pdf of a uniform lapse probability distribution on the response range
    randomPDF = 1/(RespRange(2)-RespRange(1));                              
    if (nTrialsA > 0); LLvectorA(:) = log(randomPDF*ones(nTrialsA,1)); end                   
    if (nTrialsV > 0); LLvectorV(:) = log(randomPDF*ones(nTrialsV,1)); end
end

%Common source judgments
if isnan(sA) || isnan(sV)
    %Unisensory
    likeFunCSJ = [0 1]; %[Common = 0, Indep = 1];
else
    %Bisensory
    likeFunCSJ = (1/2)*ones(1,2);
end
if (nTrialsC > 0)
    LLvectorC(:) = log(likeFunCSJ(CSJ_resp));                               %Note that the responses were already transformed to bin_idx numbers before
end
    
end %[EOF]
