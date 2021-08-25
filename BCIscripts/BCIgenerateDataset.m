function [trueLocsAV,responsesAVC,i_Conditions] = BCIgenerateDataset(P,StimLocs,nTrialsPerAVcond)
% Generate a fake dataset for a factorial design (sA,sV) that includes 
% unisensory trials. Input "P" is a cell-array of structures with parameter 
% settings; one cell per experimental condition. 

%Set some defaults
if nargin < 2
    StimLocs = [-10 -5 0 5 10];
end
if nargin < 3
    nTrialsPerAVcond = 40;
end

%Find number of experimental conditions
nExpCons = length(P);

%Set all sA and sV combinations (factorial design incl. unisensory conditions)
sAs = [StimLocs,NaN];
sVs = [StimLocs,NaN];
nAVcons = numel(sAs)*numel(sVs)-1;                              
sAV = nan(nAVcons,2);
for i=1:numel(sAs)
    for j=1:numel(sVs)
        if ~(i==numel(sAs) && j==numel(sVs))
            idxTmp = (i-1)*numel(sVs)+j;
            sAV(idxTmp,:) = [sAs(i) sVs(j)];
        end
    end
end

%Generate one dataset by looping over all sA,sV location combinations
nTrialsTotal = nExpCons*nAVcons*nTrialsPerAVcond;
i_Conditions = false(nTrialsTotal,nExpCons);
trueLocsAV = nan(nTrialsTotal,2);
responsesAVC = nan(nTrialsTotal,3);
for c=1:nExpCons
    for i=1:nAVcons
        index = (c-1)*nAVcons*nTrialsPerAVcond + (i-1)*nTrialsPerAVcond + (1:nTrialsPerAVcond);
        i_Conditions(index,c) = true;
        trueLocsAV(index,:) = repmat(sAV(i,:),[nTrialsPerAVcond 1]);
        [responsesAVC(index,1), responsesAVC(index,2), responsesAVC(index,3)] = BCIgenerateResponsesMonteCarlo(P{c},sAV(i,1),sAV(i,2),nTrialsPerAVcond);
    end
end

%Randomize the order of the trials
mixOrder = randperm(nTrialsTotal);
trueLocsAV = trueLocsAV(mixOrder,:);
responsesAVC = responsesAVC(mixOrder,:);
i_Conditions = i_Conditions(mixOrder,:);

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function [sA_resp, sV_resp, CSJ_resp] = BCIgenerateResponsesMonteCarlo(P,sA,sV,nSamples)
% Generate BCI responses for one condition (combination of sA and sV)

if nargin < 4
    nSamples = 40;
end

% Set defaults
RespLocs = NaN;                     %Response Locations: E.g. NaN or [-10 -5 0 5 10] for continuous or discrete responses respectively
RespRange = [-20 20];               %Response range in case of continuous responses (e.g. screen limits)

Pcommon = 0.5;                      %common source prior
sigmaA = 5;                         %sigma of auditory noise
sigmaV = 2;                         %sigma of visual noise
sigmaP = 10;                        %sigma of spatial prior
muP = 0;                            %centre of spatial prior

deltaSigmaA = 0;                    %Linear sigmaA shift with eccentricity: sigmaA_final = sigmaA + deltaSigmaA*sA  
deltaSigmaV = 0;                    %Linear sigmaV shift with eccentricity: sigmaV_final = sigmaV + deltaSigmaV*sV  
deltaXA = 0;                        %Linear xA shift with eccentricity: xA_final = xA + deltaXA*xA 
deltaXV = 0;                        %Linear xV shift with eccentricity: xV_final = xV + deltaXV*xV 

sigmaMotor = 1;                     %sigma of motor noise (only used in case of continuous responses)

decisionFun = 'ModelAveraging';     %BCI decision function (see Wozny, Beierholm, Shams, 2010, PLOS Comp. Biology): 'ModelAveraging', 'ModelSelection' or 'ProbabilityMatching'
CSJthresh = 0.5;                    %common source judgment threshold (only used with ModelAveraging decision function)

lapseRA = max(1/nSamples,0.02);     %Aud
lapseRV = max(1/nSamples,0.02);     %Vis
lapseRC = max(1/nSamples,0.02);     %CSJ
lapseR = NaN;                       %lapse rate for all responses (i.e. Overwrites others if set to non-NaN)

% Overwrite defaults with entries from struct "P"
% Let OptionsStruct overwrite the default values (defined by the user)
pFields = fieldnames(P);
for i=1:length(pFields)
    eval([pFields{i,1} ' = P.(pFields{i,1});']);    
end

% Deal with a shared lapse rate across auditory and visual modalities...
% Overwrite the modality-specific lapse rates with the shared lapse rate
if ~isnan(lapseR)
    lapseRA = lapseR;
    lapseRV = lapseR;
    lapseRC = lapseR;
end

% Correct sensory noise for eccentricity
sigmaA = max(sigmaA + deltaSigmaA*abs(sA),1e-3);
sigmaV = max(sigmaV + deltaSigmaV*abs(sV),1e-3);

% Variances of A and V and prior                                            
varA = sigmaA^2;
varV = sigmaV^2;
varP = sigmaP^2;

% Variances of posteriors: P(S|xA,xV,C=1) and P(S|xV,C=2) and P(S|xA,C=2)
varA_pos = 1/(1/varA + 1/varP);
varV_pos = 1/(1/varV + 1/varP);
varAV_pos = 1/(1/varV + 1/varA + 1/varP);

% Variances of common/segregated source(s) likelihood: P(xA,xV|C=1) and P(xA,xV|C=2) 
varA_PxGivenIndep = varA + varP;
varV_PxGivenIndep = varV + varP;
varAV_PxGivenCommon = varV*varA + varV*varP + varA*varP;

% Simulate internal representations
xV = sV + sigmaV*randn(nSamples,1);  
xA = sA + sigmaA*randn(nSamples,1); 

% Correct the likelihood distributions for their eccentricity-dependent bias
xV = xV + deltaXV*xV;
xA = xA + deltaXA*xA;

% Maximum a posteriori (MAP) estimates given common and independent causes
sV_hat_indep = (xV/varV + muP/varP) * varV_pos;
sA_hat_indep = (xA/varA + muP/varP) * varA_pos;

%Unisensory?
if isnan(sA) || isnan(sV)

    sV_hat = sV_hat_indep;
    sA_hat = sA_hat_indep;
    CSJ_resp = false(nSamples,1);
    
%Bisensory?    
else
    
    sAV_hat_common = (xV/varV + xA/varA + muP/varP) * varAV_pos;

    % The likelihoods of common source and independent sources: P(xA,xV|C=1) and P(xA,xV|C=2)
    likelihood_common = exp(-((xV-xA).^2 * varP + (xV-muP).^2 * varA + (xA-muP).^2 * varV)/(2*varAV_PxGivenCommon))/(2*pi*sqrt(varAV_PxGivenCommon));
    likelihoodV_indep = exp(-((xV-muP).^2)/(2*varV_PxGivenIndep))/sqrt(2*pi*varV_PxGivenIndep);
    likelihoodA_indep = exp(-((xA-muP).^2)/(2*varA_PxGivenIndep))/sqrt(2*pi*varA_PxGivenIndep);
    likelihood_indep =  likelihoodV_indep .* likelihoodA_indep;

    % Posterior probability of common source: P(C=1|xA,xV)
    post_common = likelihood_common * Pcommon;
    post_indep = likelihood_indep * (1-Pcommon);
    pC = post_common./(post_common + post_indep);

    % Final location estimate depends on the decision function
    if strcmpi(decisionFun,'ModelAveraging')

        % Weighted average of integration and segregation
        sV_hat = pC .* sAV_hat_common + (1-pC) .* sV_hat_indep;
        sA_hat = pC .* sAV_hat_common + (1-pC) .* sA_hat_indep;
        CSJ_resp = pC >= CSJthresh;
        
    elseif strcmpi(decisionFun,'ModelSelection')

        % Choose integration of segregation based on highest probability
        sV_hat = (pC>0.5).*sAV_hat_common + (pC<=0.5).*sV_hat_indep;
        sA_hat = (pC>0.5).*sAV_hat_common + (pC<=0.5).*sA_hat_indep;
        CSJ_resp = pC >= 0.5;
        
    elseif strcmpi(decisionFun,'ProbabilityMatching')

        % Probability matching (stochastically choose integration or segregation, where the probability for either choice is defined by P(C=1|xA,xV)  
        alpha = rand(nSamples,1); 
        sV_hat = (pC>alpha).*sAV_hat_common + (pC<=alpha).*sV_hat_indep;
        sA_hat = (pC>alpha).*sAV_hat_common + (pC<=alpha).*sA_hat_indep;
        CSJ_resp = pC >= alpha;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Discretize for response buttons (?) and add lapse rate %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Select at random a certain number of generated samples that are lapses. We overwrite those simulated responses..    
nLapseSamplesV = floor(lapseRV*nSamples) + (mod(lapseRV*nSamples,1) >= rand(1));
nLapseSamplesA = floor(lapseRA*nSamples) + (mod(lapseRA*nSamples,1) >= rand(1));
nLapseSamplesC = floor(lapseRC*nSamples) + (mod(lapseRC*nSamples,1) >= rand(1));
idx_lapseV = datasample(1:nSamples,nLapseSamplesV,'Replace',false);
idx_lapseA = datasample(1:nSamples,nLapseSamplesA,'Replace',false);
idx_lapseC = datasample(1:nSamples,nLapseSamplesC,'Replace',false);

% Discretize for response buttons?
if all(~isnan(RespLocs))
    
    %Find the response button closest to sV_hat and sA_hat and choose those response buttons    
    nRespLocs = numel(RespLocs);
    
    if ~isnan(sV)
        [~,idxV] = min(abs(repmat(sV_hat,[1,nRespLocs]) - repmat(RespLocs,[nSamples,1])),[],2);
        sV_resp = RespLocs(idxV)';
        sV_resp(idx_lapseV,1) = RespLocs(randi(nRespLocs,[nLapseSamplesV,1]));      %Lapses overwrite some (randomly chosen) other responses
    else
        sV_resp = nan(nSamples,1);                                                  %Generate no visual responses if there was no visual stimulus
    end
    
    if ~isnan(sA)
        [~,idxA] = min(abs(repmat(sA_hat,[1,nRespLocs]) - repmat(RespLocs,[nSamples,1])),[],2);
        sA_resp = RespLocs(idxA)';
        sA_resp(idx_lapseA,1) = RespLocs(randi(nRespLocs,[nLapseSamplesA,1]));      %Lapses overwrite some (randomly chosen) other responses
    else
        sA_resp = nan(nSamples,1);                                                  %Generate no auditory responses if there was no auditory stimulus
    end
    
%Don't discretize if RespLocs is NaN (continuous responses) 
%Lapses are chosen from a uniform distribution on a certain response range, such that the distribution is completely independent of the other parameters in the BCI model
else
    %Add motor noise to the continuous responses
    sV_hat = sV_hat + sigmaMotor*randn(nSamples,1);
    sA_hat = sA_hat + sigmaMotor*randn(nSamples,1);
    
    %Ensure that responses are not out of response range
    sV_hat(sV_hat < RespRange(1)) = RespRange(1);
    sV_hat(sV_hat > RespRange(2)) = RespRange(2);
    sA_hat(sA_hat < RespRange(1)) = RespRange(1);
    sA_hat(sA_hat > RespRange(2)) = RespRange(2);
    
    sV_resp = sV_hat;        %No need for overwriting unisensory conditions with NaNs, because sV_hat and sA_hat are already NaNs if respective sV or sA are NaN
    if ~isnan(sV)
        sV_resp(idx_lapseV,1) = (RespRange(2)-RespRange(1))*rand([nLapseSamplesV,1])+RespRange(1);      %Lapses overwrite some (randomly chosen) other responses
    end
    
    sA_resp = sA_hat;
    if ~isnan(sA)
        sA_resp(idx_lapseA,1) = (RespRange(2)-RespRange(1))*rand([nLapseSamplesA,1])+RespRange(1);      %Lapses overwrite some (randomly chosen) other responses
    end
    
end

%Overwrite some common source judgment with random lapses (if not unisensory)
if ~(isnan(sA) || isnan(sV))
    PossibleResponses = [false true];
    CSJ_resp(idx_lapseC,1) = PossibleResponses(randi(2,[nLapseSamplesC,1]));
end

end %[EOF]
