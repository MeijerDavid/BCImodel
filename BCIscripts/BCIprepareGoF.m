function [LLvectors,InfoVectors,countsResp,PredResps_allCons] = BCIprepareGoF(paramsConv,P,trueLocsAV,responsesAVC,i_Conditions)
%Compute the log-likelihoods of the BCI model, random (null) model and the maximum model, return them as column vectors 
%that contain the LL of each response (with indices indicating whether they are auditory, visual, or common-source responses).    
%Count the number of discrete responses per RespLoc and per location/condition.   
%In addition, save all BCI predicted responses (in a cell array, per condition, and then per AV location combination).

%%% - The build-up of this script is very similar to BCIcompLL.m - %%%

%Convert the parameters back to 'normal' space
paramsReal = BCIconvertParams2Fit(paramsConv,P,'conv2real');

%Gather all unique locations (combinations of sA and sV)
LocConds = P.sAsVconds;
nLocConds = size(LocConds,1);

%Initialize matrices wherein we save the LL of each response for the three models: BCI, random and max   
LLresponsesBCI = nan(size(responsesAVC));
LLresponsesRandom = nan(size(responsesAVC));
LLresponsesMax = nan(size(responsesAVC));

%Keep track of other conditions too (relating to the LLs above)
ConIdx = nan(size(responsesAVC));
Alocs  = nan(size(responsesAVC));
Vlocs  = nan(size(responsesAVC));
Aresps = nan(size(responsesAVC));
Vresps = nan(size(responsesAVC));
Cresps = nan(size(responsesAVC));

%Initialize a cell array to store all predicted response distributions etc.
PredResps_allCons = cell(nLocConds,P.nConditions);

%Initialize response counts for discrete responses
responseCountsCSJ = zeros(P.nConditions*nLocConds,2);
responseCountsLocA = zeros(P.nConditions*nLocConds,numel(P.RespLocs));
responseCountsLocV = zeros(P.nConditions*nLocConds,numel(P.RespLocs));
responseCountsInfo.ConIdx = nan(P.nConditions*nLocConds,1);
responseCountsInfo.Alocs = nan(P.nConditions*nLocConds,1);
responseCountsInfo.Vlocs = nan(P.nConditions*nLocConds,1);
if all(~isnan(P.RespLocs)) %If discrete responses
    discreteBool = 1;
    %create bins to use for the histogram (we use 'histcounts' to avoid problems with numerical approximations)
    if numel(P.RespLocs) == 1
        edges = [-inf inf]; 
    else
        edges = [-inf mean([P.RespLocs(1:(end-1)); P.RespLocs(2:end)],1) inf];  
    end
else %Continuous responses
    discreteBool = 0;
end

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
        PredResps_allCons{i,c}.sA = sA;
        PredResps_allCons{i,c}.sV = sV;
        
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
            sA_resp = responsesAVC(i_relResp_A,1);                                              %These can be 'empty'
            sV_resp = responsesAVC(i_relResp_V,2);
            CSJ_resp = responsesAVC(i_relResp_C,3);
            
            %If there are valid responses...
            if sum([numel(sA_resp); numel(sV_resp); numel(CSJ_resp)]) > 0   
                
                %Save the actual responses for this condition    
                PredResps_allCons{i,c}.respA = responsesAVC(i_AVcond,1);                        %These may contain NaNs and will have been indexed for discrete responses
                PredResps_allCons{i,c}.respV = responsesAVC(i_AVcond,2);                               
                PredResps_allCons{i,c}.respCSJ = responsesAVC(i_AVcond,3);
                
                %Compute log-likelihoods of the responses for the three models: BCI, random and maxLL   
                [LLvectorA_bci, LLvectorV_bci, LLvectorC_bci, likeFunA, likeFunV, likeFunCSJ, likeFunA_grid, likeFunV_grid] = BCIcompLL_bci(Pcopy,sA,sV,sA_resp,sV_resp,CSJ_resp);
                [LLvectorA_random, LLvectorV_random, LLvectorC_random] = BCIcompLL_random(Pcopy,sA,sV,sA_resp,sV_resp,CSJ_resp);
                [LLvectorA_max, LLvectorV_max, LLvectorC_max] = BCIcompLL_max(Pcopy,sA,sV,sA_resp,sV_resp,CSJ_resp);
            
                %Save the predicted response frequencies of the BCI model as a structure   
                PredResps_allCons{i,c}.likeFunA = likeFunA;
                PredResps_allCons{i,c}.likeFunV = likeFunV;
                PredResps_allCons{i,c}.likeFunCSJ = likeFunCSJ;
                PredResps_allCons{i,c}.likeFunA_grid = likeFunA_grid;
                PredResps_allCons{i,c}.likeFunV_grid = likeFunV_grid;
            
                %Save the log-likelihood of each response in the matrix
                LLresponsesBCI(i_relResp_A,1) = LLvectorA_bci;          LLresponsesBCI(i_relResp_V,2) = LLvectorV_bci;          LLresponsesBCI(i_relResp_C,3) = LLvectorC_bci;
                LLresponsesRandom(i_relResp_A,1) = LLvectorA_random;    LLresponsesRandom(i_relResp_V,2) = LLvectorV_random;    LLresponsesRandom(i_relResp_C,3) = LLvectorC_random;            
                LLresponsesMax(i_relResp_A,1) = LLvectorA_max;          LLresponsesMax(i_relResp_V,2) = LLvectorV_max;          LLresponsesMax(i_relResp_C,3) = LLvectorC_max;

                %Save the other conditional information 
                ConIdx(i_relResp_A,1) = c;                              ConIdx(i_relResp_V,2) = c;                              ConIdx(i_relResp_C,3) = c;
                Alocs(i_relResp_A,1) = sA;                              Alocs(i_relResp_V,2) = sA;                              Alocs(i_relResp_C,3) = sA;
                Vlocs(i_relResp_A,1) = sV;                              Vlocs(i_relResp_V,2) = sV;                              Vlocs(i_relResp_C,3) = sV;   
                if discreteBool
                    Aresps(i_relResp_A,1) = P.RespLocs(sA_resp);        Vresps(i_relResp_V,2) = P.RespLocs(sV_resp);            Cresps(i_relResp_C,3) = CSJ_resp;       
                else
                    Aresps(i_relResp_A,1) = sA_resp;                    Vresps(i_relResp_V,2) = sV_resp;                        Cresps(i_relResp_C,3) = CSJ_resp;       
                end
            end
            
            %Count the participants responses (to store in the countsArray, see further below)    
            rowIndex = (c-1)*nLocConds+i;
            responseCountsInfo.ConIdx(rowIndex,1) = c;
            responseCountsInfo.Alocs(rowIndex,1) = sA;
            responseCountsInfo.Vlocs(rowIndex,1) = sV;
            if discreteBool
                responseCountsLocA(rowIndex,:) = histcounts(P.RespLocs(sA_resp),edges);
                responseCountsLocV(rowIndex,:) = histcounts(P.RespLocs(sV_resp),edges);
            end
            n1 = sum(CSJ_resp == 1);
            n2 = sum((CSJ_resp ~= 1) & ~isnan(CSJ_resp));
            responseCountsCSJ(rowIndex,:) = [n1 n2];
            
        end
    end %end for-loop of locations
    
end %end for-loop of conditions

%Get the relevant indices of the responses
i_rel1 = ~isnan(LLresponsesBCI(:,1));         i_rel2 = ~isnan(LLresponsesBCI(:,2));         i_rel3 = ~isnan(LLresponsesBCI(:,3)); 

%Create long vectors with log-likelihoods for all responses (ignoring NaNs)   
LLvectors.bci     = [LLresponsesBCI(i_rel1,1);        LLresponsesBCI(i_rel2,2);         LLresponsesBCI(i_rel3,3)];
LLvectors.random  = [LLresponsesRandom(i_rel1,1);     LLresponsesRandom(i_rel2,2);      LLresponsesRandom(i_rel3,3)];
LLvectors.max     = [LLresponsesMax(i_rel1,1);        LLresponsesMax(i_rel2,2);         LLresponsesMax(i_rel3,3)];

%Do the same for the other conditional info
InfoVectors.ConIdx    = [ConIdx(i_rel1,1);        ConIdx(i_rel2,2);           ConIdx(i_rel3,3)];
InfoVectors.Alocs     = [Alocs(i_rel1,1);         Alocs(i_rel2,2);            Alocs(i_rel3,3)];
InfoVectors.Vlocs     = [Vlocs(i_rel1,1);         Vlocs(i_rel2,2);            Vlocs(i_rel3,3)];
InfoVectors.Aresps    = [Aresps(i_rel1,1);        Aresps(i_rel2,2);           Aresps(i_rel3,3)];
InfoVectors.Vresps    = [Vresps(i_rel1,1);        Vresps(i_rel2,2);           Vresps(i_rel3,3)];
InfoVectors.CSJresps  = [Cresps(i_rel1,1);        Cresps(i_rel2,2);           Cresps(i_rel3,3)];

%Create the response counts arrays...
i_noAresps = (sum(responseCountsLocA,2) == 0);
countsResp.A.counts = responseCountsLocA(~i_noAresps,:);
countsResp.A.info.ConIdx = responseCountsInfo.ConIdx(~i_noAresps,1);
countsResp.A.info.Alocs = responseCountsInfo.Alocs(~i_noAresps,1);
countsResp.A.info.Vlocs = responseCountsInfo.Vlocs(~i_noAresps,1);

i_noVresps = (sum(responseCountsLocV,2) == 0);
countsResp.V.counts = responseCountsLocV(~i_noVresps,:);
countsResp.V.info.ConIdx = responseCountsInfo.ConIdx(~i_noVresps,1);
countsResp.V.info.Alocs = responseCountsInfo.Alocs(~i_noVresps,1);
countsResp.V.info.Vlocs = responseCountsInfo.Vlocs(~i_noVresps,1);

i_noCSJresps = (sum(responseCountsCSJ,2) == 0);
countsResp.CSJ.counts = responseCountsCSJ(~i_noCSJresps,:);
countsResp.CSJ.info.ConIdx = responseCountsInfo.ConIdx(~i_noCSJresps,1);
countsResp.CSJ.info.Alocs = responseCountsInfo.Alocs(~i_noCSJresps,1);
countsResp.CSJ.info.Vlocs = responseCountsInfo.Vlocs(~i_noCSJresps,1);

end %[EOF]
