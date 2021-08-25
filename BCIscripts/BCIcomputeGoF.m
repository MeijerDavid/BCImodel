function GoF = BCIcomputeGoF(P, trueLocsAV, responsesAVC, i_Conditions, BadsOrMcmc, varargin)
%Compute the goodness-of-fit measures

%BADS
if strcmp(BadsOrMcmc,'BADS') 
    
    fittedParamsBADS = varargin{1};
    
    %Compute log-likelihood vectors (one entry for each response) of the BCI model, random (null) model and the maximum model
    %Also return vectors with associated conditional and response info.    
    %Count the number of discrete responses per RespLoc and per location/condition.   
    %And save all BCI predicted responses (in a cell array, per condition, and then per AV location combination).
    [LLvectors,InfoVectors,countsResp,PredResps_allCons] = BCIprepareGoF(fittedParamsBADS,P,trueLocsAV,responsesAVC,i_Conditions);
    
    %Set up output structure
    GoF = [];
    GoF.discreteCounts = countsResp;
    GoF.predictedResps = PredResps_allCons;
    GoF.LLvectors = LLvectors;
    GoF.InfoVectors = InfoVectors;

%MCMC (compute R2 and absGoF again using the cross-validated LLs)   
elseif strcmp(BadsOrMcmc,'MCMC') 
    
    %Retrieve the old GoF info from the input
    GoF_bads = varargin{1};
    mcmcLL_responses = varargin{2};
    minNeff = varargin{3};
    
    %Compute relative MCMC efficiency N_eff/N 
    Reff = minNeff/size(mcmcLL_responses,1);
    
    %There's no need to compute the minimum and maximum LL again! (Just overwrite the LLs for the BCI model with the cross-validated LLs from MCMC, see below)   
    LLvectors = GoF_bads.LLvectors;
    InfoVectors = GoF_bads.InfoVectors;
    countsResp = GoF_bads.discreteCounts;
    
    %Call Aki Vehtari's psisloo function (Pareto smoothed importance sampling leave-one-out log predictive densities) 
    %https://github.com/avehtari/PSIS/blob/master/m/psisloo.m
    [loo,loos,pk] = psisloo_DM(mcmcLL_responses,Reff);
    PSISLOO.loo = loo;                              %The leave-one-out log predictive density 
    PSISLOO.loos = loos;                            %The individual leave-one-out log predictive densities (per data point)
    PSISLOO.kIdx = pk;                              %The Pareto tail k indices
    PSISLOO.kn0 = sum(pk<0.5);                      %The number of data points for which the estimation was successful
    PSISLOO.kn1 = sum(pk>=0.5&pk<1);                %These numbers should be small!
    PSISLOO.kn2 = sum(pk>=1);                       %These numbers should be small!  
    
    %Correct the BCI log-likelihoods   
    LLvectors.bci = PSISLOO.loos;
    
    %Set up output structure
    GoF = [];
    GoF.PSISLOO = PSISLOO;   
end

%Discrete localisation responses?
if all(~isnan(P.RespLocs)) 
    discreteBool = 1;
else
    discreteBool = 0;
end

%Initialize R2 and absGoF cell arrays
if P.nConditions > 1
    
    ConNames = cell(1,P.nConditions+1);
    ConNames{1,1} = 'Overall';
    iCon = cell(4,P.nConditions+1);
    iCon{1,1} = true(numel(LLvectors.bci),1);                   %LLs
    iCon{2,1} = true(numel(countsResp.A.info.ConIdx),1);        %counts.A
    iCon{3,1} = true(numel(countsResp.V.info.ConIdx),1);        %counts.V
    iCon{4,1} = true(numel(countsResp.CSJ.info.ConIdx),1);      %counts.CSJ
    for c=1:P.nConditions 
        ConNames{1,c+1} = ['Con ' num2str(c)]; 
        iCon{1,c+1} = InfoVectors.ConIdx == c;                  %LLs
        iCon{2,c+1} = countsResp.A.info.ConIdx == c;            %counts.A
        iCon{3,c+1} = countsResp.V.info.ConIdx == c;            %counts.V
        iCon{4,c+1} = countsResp.CSJ.info.ConIdx == c;          %counts.CSJ
    end
    
    R2 = cell(2,P.nConditions+1);
    R2(1,:) = ConNames;
    absGoF = cell(2,P.nConditions+1);
    absGoF(1,:) = ConNames;
    
else %only one "condition"
    iCon = {true(numel(LLvectors.bci),1); true(numel(countsResp.A.info.ConIdx),1); true(numel(countsResp.V.info.ConIdx),1); true(numel(countsResp.CSJ.info.ConIdx),1)};
    R2 = cell(2,1);
    R2{1,1} = 'Overall';
    absGoF = cell(2,1);
    absGoF{1,1} = 'Overall';
end
    
%Compute the goodness-of-fit measures (R2 and absGoF) per condition
for i=1:size(iCon,2)
    
    %Initialize
    absGoFtmp = [];
    R2tmp = [];
    
    %Gather the positional info
    i_total   = iCon{1,i};
    i_respCSJ_only = iCon{1,i} & ~isnan(InfoVectors.CSJresps);
    i_respA_only   = iCon{1,i} & ~isnan(InfoVectors.Aresps);
    i_respV_only   = iCon{1,i} & ~isnan(InfoVectors.Vresps);
    i_respA_and_respV  = i_respA_only | i_respV_only; 
    
    %Compute the total R-squared GoF measure based on the method by Nagelkerke, 1991, Biometrika (see helper function below)    
    R2tmp.total = computeR2(sum(LLvectors.bci(i_total)),sum(LLvectors.random(i_total)),sum(LLvectors.max(i_total)),sum(i_total));
    
    %All localisation responses (A and V, but not CSJ)
    if sum(i_respA_and_respV) > 0
        R2tmp.respA_and_respV = computeR2(sum(LLvectors.bci(i_respA_and_respV)),sum(LLvectors.random(i_respA_and_respV)),sum(LLvectors.max(i_respA_and_respV)),sum(i_respA_and_respV));
        %For discrete responses...
        if discreteBool 
            %Compute absolute GoF based on the Kullback-Leibler Divergence (Acerbi, Dokka, Angelaki & Ma, 2018, PLOS comp. Biology). See: https://github.com/lacerbi/gofit
            countsMatrix = [countsResp.A.counts(iCon{2,i},:); countsResp.V.counts(iCon{3,i},:)];
            absGoFtmp.respA_and_respV = gofit(countsMatrix,sum(LLvectors.bci(i_respA_and_respV)));
        end
    end
    
    %Auditory localisation responses only (A)
    if sum(i_respA_only) > 0
        R2tmp.respA_only = computeR2(sum(LLvectors.bci(i_respA_only)),sum(LLvectors.random(i_respA_only)),sum(LLvectors.max(i_respA_only)),sum(i_respA_only));
        if discreteBool 
            absGoFtmp.respA_only = gofit(countsResp.A.counts(iCon{2,i},:),sum(LLvectors.bci(i_respA_only)));
        end
    end
    
    %Visual localisation responses only (V)
    if sum(i_respV_only) > 0
        R2tmp.respV_only = computeR2(sum(LLvectors.bci(i_respV_only)),sum(LLvectors.random(i_respV_only)),sum(LLvectors.max(i_respV_only)),sum(i_respV_only));
        if discreteBool 
            absGoFtmp.respV_only = gofit(countsResp.V.counts(iCon{3,i},:),sum(LLvectors.bci(i_respV_only)));
        end
    end
    
    %Common source judgments only (CSJ)
    if sum(i_respCSJ_only) > 0
        R2tmp.respCSJ_only = computeR2(sum(LLvectors.bci(i_respCSJ_only)),sum(LLvectors.random(i_respCSJ_only)),sum(LLvectors.max(i_respCSJ_only)),sum(i_respCSJ_only));
        absGoFtmp.respCSJ_only = gofit(countsResp.CSJ.counts(iCon{4,i},:),sum(LLvectors.bci(i_respCSJ_only)));
    end
    
    %Save outputs for this condition
    R2{2,i} = R2tmp;
    absGoF{2,i} = absGoFtmp;

end %end for-loop across conditions

%Gather remainder of output structure
GoF.R2 = R2;
GoF.absGoF = absGoF;

end %[EOF]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function %%%
%%%%%%%%%%%%%%%%%%%%%%%

%Compute R-squared measures (based on likelihood ratios)
function R2 = computeR2(LLbci,LLrandom,LLmax,n)
    BCIvsNull = 1 - exp(-(2/n)*(LLbci-LLrandom));                           %Compare BCI vs Null model (random responses). This R-squared can be negative if the BCI model performs worse than the random model!
    max_abs = 1 - exp((2/n)*LLrandom);                                      %The absolute maximally achievable R2 (if all responses were deterministically predicted with probability of 1; i.e. LLbci = 0)
    max_prob = 1 - exp(-(2/n)*(LLmax-LLrandom));                            %The maximally achievable R2 for a probabilistic/stochastic/non-deterministic model such as the BCI model
    R2.BCIvsNull_Nagelkerke = BCIvsNull/max_abs;                            %R-squared normalized relative to the absolute maximum R2, as proposed by Nagelkerke (1991, Biometrika)                           
    R2.BCIvsNull_Stochastic = BCIvsNull/max_prob;                           %R-squared normalized relative to the maximum R2 for the current type of model (i.e. R2 = 1 here truly means that the goodness-of-fit is maximal)
    
    R2.helper.BCIvsNull = BCIvsNull;                                        %Save the 'helper R2' values (see explanations above)
    R2.helper.max_abs = max_abs;
    R2.helper.max_prob = max_prob;
end
