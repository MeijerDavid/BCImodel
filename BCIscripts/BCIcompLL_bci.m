function [LLvectorA,LLvectorV,LLvectorC,likeFunA,likeFunV,likeFunCSJ,likeFunA_grid,likeFunV_grid] = BCIcompLL_bci(P,sA,sV,sA_resp,sV_resp,CSJ_resp,randomNoiseFlag)
% Generate BCI response frequencies for one condition (combination of sA and sV)
% And compute loglike for A, V and CSJ responses

if nargin < 7
    randomNoiseFlag = 0;    %Default: no random noise added to numerical integration method.
end                         %With random noise is theoretically better, but causes more difficulties for BADS, VBMC and MCMC algorithms.

% Gather relevant variables from structure P
nGrid = P.nGridNumInt;
triangpdfSpeedUp = P.triangpdfSpeedUp;
integrateMethod = P.integrateMethod;

RespLocs = P.RespLocs;
RespRange = P.RespRange;

Pcommon = P.Pcommon;
sigmaA = P.sigmaA;
sigmaV = P.sigmaV;
sigmaP = P.sigmaP;
muP = P.muP;

deltaSigmaA = P.deltaSigmaA;                    
deltaSigmaV = P.deltaSigmaV;                    
deltaXA = P.deltaXA;                        
deltaXV = P.deltaXV;                        

sigmaMotor = P.sigmaMotor;
decisionFun = P.decisionFun;
CSJthresh = P.CSJthresh;

lapseR = P.lapseR;
lapseRA = P.lapseRA;
lapseRV = P.lapseRV;
lapseRC = P.lapseRC;

% Deal with a shared lapse rate across auditory and visual modalities...
% Overwrite the modality-specific lapse rates with the shared lapse rate
if ~isnan(lapseR)
    lapseRA = lapseR;
    lapseRV = lapseR;
    lapseRC = lapseR;
end

% Determine number of trials
nTrialsA = numel(sA_resp);
nTrialsV = numel(sV_resp);
nTrialsC = numel(CSJ_resp);

% Initialize output vectors
LLvectorA = nan(nTrialsA,1);
LLvectorV = nan(nTrialsV,1);
LLvectorC = nan(nTrialsC,1);

% Correct sensory noise for eccentricity
sigmaA = sigmaA + deltaSigmaA*abs(sA-muP);
sigmaV = sigmaV + deltaSigmaV*abs(sV-muP);

% Prepare the grids for numerical integration
gridCoverSD = 6;                                                            %Ensure that we integrate over 'all' xA / xV (from -3*SD to +3*SD)
A_spacing = gridCoverSD*sigmaA/(nGrid-1);
V_spacing = gridCoverSD*sigmaV/(nGrid-1);
noiseShiftA = randomNoiseFlag*(rand(1)-0.5)*A_spacing;                      %Add random noise to centering of grid
noiseShiftV = randomNoiseFlag*(rand(1)-0.5)*V_spacing;

% Initialize some settings for discrete responses
if all(~isnan(RespLocs))
    discreteBool = 1;
    nRespLocs = numel(RespLocs);
    likeFunA = nan(1,nRespLocs);                                            %Initialize (may not change if unisensory condition or no response in this modality)
    likeFunV = nan(1,nRespLocs);    
    likeFunA_grid = RespLocs;
    likeFunV_grid = RespLocs;
    if nRespLocs == 1
        Resp_edges_high = inf;
    else
        Resp_edges_mean = (RespLocs(1:(end-1)) + RespLocs(2:end))/2;
        Resp_edges_high = [Resp_edges_mean inf];
    end
% Initialize some settings for continuous responses   
else
    discreteBool = 0;
    likeFunA = NaN; likeFunA_grid = NaN;                                    %Initialize (may not change if unisensory condition or no response in this modality)
    likeFunV = NaN; likeFunV_grid = NaN;  
    lapsePDF = 1/(RespRange(2)-RespRange(1));                               %The constant pdf of a uniform lapse probability distribution on the response range
    if integrateMethod == 1
        if triangpdfSpeedUp                                                 %To use, set: OptionsStruct.triangpdfSpeedUp = 1;
            motorNoisePDFfun = @(x,mu,sigma) bsxfun_triangpdf(x,mu,sigma);  %Triangular motor noise pdf as an approximation to gaussian motor noise pdf (about 3 times faster than "bsxfun_normpdf")  
        else
            motorNoisePDFfun = @(x,mu,sigma) bsxfun_normpdf(x,mu,sigma);    %Default Gaussian motor noise pdf
        end
    end
end
likeFunCSJ = nan(1,2);                                                      %Initialize common source judgment likelihood function

% Pre-compute some variances   
varA = sigmaA^2;                                                            %Variances of A and V and prior                                            
varV = sigmaV^2;
varP = sigmaP^2;

varA_pos = 1/(1/varA + 1/varP);                                             %Variances of posteriors: P(S|xA,xV,C=1) and P(S|xV,C=2) and P(S|xA,C=2)
varV_pos = 1/(1/varV + 1/varP);
varAV_pos = 1/(1/varV + 1/varA + 1/varP);

varA_PxGivenIndep = varA + varP;                                            %Variances of common/segregated source(s) likelihood: P(xA,xV|C=1) and P(xA,xV|C=2) 
varV_PxGivenIndep = varV + varP;
varAV_PxGivenCommon = varV*varA + varV*varP + varA*varP;

% Unisensory auditory stimulus (A)
if isnan(sV)
    
    % Auditory responses
    if nTrialsA > 0
    
        % Sample xA on regular grid and compute the approximate probabilities of the bins   
        xA = linspace(sA-(gridCoverSD/2)*sigmaA,sA+(gridCoverSD/2)*sigmaA,nGrid)' + noiseShiftA;        %Sample xA on a regular grid (column vector)
        xA_prob = normpdf(xA,sA,sigmaA)*A_spacing;                                                      %Compute the probabilities for each xA    

        % Correct the likelihood distribution for its eccentricity-dependent bias
        xA = xA + deltaXA*(xA-muP);

        % Bias by prior: compute unisensory maximum a posteriori (MAP) estimates for each xA   
        sA_hat = (xA/varA + muP/varP) * varA_pos;

        % Discrete responses
        if discreteBool

            %Use simple trapezoidal numerical integration
            cdf_likeA = qtrapz(bsxfun(@times,xA_prob,bsxfun(@ge,Resp_edges_high,sA_hat)),1);            %Integrate out xA (1st dim)   
            likeFunA = diff([0 reshape(cdf_likeA,[1 nRespLocs])]);                                      %Convert back to probabilities  
            likeFunA = lapseRA*(1/nRespLocs) + (1-lapseRA)*likeFunA;                                    %Correct the likelihood function for the lapse rate
            LLvectorA(:) = log(likeFunA(sA_resp));                                                      %Save the log-likelihoods per response                             
                                                                                                        %Note that sA_resp was already transformed to bin_idx numbers before    
        % Default integration method for continuous responses     
        elseif integrateMethod == 1

            % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.
            % For each "response bin" (with centres defined by sA_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sA_hat
            % We can then compute the overall likelihood for each response by integrating out the xAs (1st dim) (i.e. the likelihood is a weighted average pdf of the motor noise Gaussian) 
            likeA = qtrapz(bsxfun(@times,xA_prob,motorNoisePDFfun(sA_resp', sA_hat, sigmaMotor)),1);
            LLvectorA(:) = log(lapseRA*lapsePDF + (1-lapseRA)*likeA);                                   %Correct the likelihoods for the lapse rate pdf, log transform and save

            % If requested as output: construct an approximate likelihood function following the alternative integration method (2. see below)     
            if nargout > 3
                [likeFunA,likeFunA_grid] = BCIcompLikeFun(sA_hat,xA_prob,nGrid,RespRange,1,sigmaMotor);
                likeFunA = lapseRA*lapsePDF + (1-lapseRA)*likeFunA;                                     %Correct the likelihood function for the lapse rate
            end

        % Alternative integration method for continuous responses    
        elseif integrateMethod == 2

            % Integrate out xA using a custom-made numerical intergration function that takes into account the width of each sA_hat bin in the xA direction
            [likeFunA,likeFunA_grid,likeFunA_spacing] = BCIcompLikeFun(sA_hat,xA_prob,nGrid,RespRange,1,sigmaMotor);
            likeFunA = lapseRA*lapsePDF + (1-lapseRA)*likeFunA;                                         %Correct the likelihood function for the lapse rate
            LLvectorA(:) = log(lininterp1(likeFunA_grid,likeFunA,sA_resp,[],likeFunA_spacing));         %Interpolate likelihood function to find pdf value for each response 
        end
    end
    
    % Visual responses (for unisensory auditory stimuli) --> Set likelihood equal to lapse rate   
    if nTrialsV > 0 
        if discreteBool
            LLvectorV(:) = log(1/nRespLocs);
        else
            LLvectorV(:) = log(lapsePDF);
        end
    end                                                                                             
    
    % Common source judgments (for unisensory auditory stimuli)  
    if nTrialsC > 0 
        likeFunCSJ = lapseRC*0.5 + (1-lapseRC)*[0 1];                                               %[Common = 0, Indep = 1];
        LLvectorC(:) = log(likeFunCSJ(CSJ_resp));                                                   %Save the log-likelihoods per response                             
    end                                                                                             %Note that CSJ_resp was already transformed to bin_idx numbers before
    
% Unisensory visual stimulus (V)
elseif isnan(sA) 
    
    % Visual responses
    if nTrialsV > 0
    
        % Sample xV on regular grid and compute the approximate probabilities of the bins   
        xV = linspace(sV-(gridCoverSD/2)*sigmaV,sV+(gridCoverSD/2)*sigmaV,nGrid) + noiseShiftV;         %Sample xV on a regular grid (row vector)
        xV_prob = normpdf(xV,sV,sigmaV)*V_spacing;                                                      %Compute the probabilities for each xV  

        % Correct the likelihood distribution for its eccentricity-dependent bias
        xV = xV + deltaXV*(xV-muP);

        % Bias by prior: compute unisensory maximum a posteriori (MAP) estimates for each xV   
        sV_hat = (xV/varV + muP/varP) * varV_pos;

        % Discrete responses
        if discreteBool

            %Use simple trapezoidal numerical integration
            cdf_likeV = qtrapz(bsxfun(@times,xV_prob,bsxfun(@ge,Resp_edges_high',sV_hat)),2);           %Integrate out xV (2nd dim) 
            likeFunV = diff([0 reshape(cdf_likeV,[1 nRespLocs])]);                                      %Convert back to probabilities  
            likeFunV = lapseRV*(1/nRespLocs) + (1-lapseRV)*likeFunV;                                    %Correct the likelihood function for the lapse rate
            LLvectorV(:) = log(likeFunV(sV_resp));                                                      %Save the log-likelihoods per response                             
                                                                                                        %Note that sV_resp was already transformed to bin_idx numbers before    
        % Default integration method for continuous responses     
        elseif integrateMethod == 1

            % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.
            % For each "response bin" (with centres defined by sV_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sV_hat
            % We can then compute the overall likelihood for each response by integrating out the xVs (2nd dim) (i.e. the likelihood is a weighted average pdf of the motor noise Gaussian) 
            likeV = qtrapz(bsxfun(@times,xV_prob,motorNoisePDFfun(sV_resp, sV_hat, sigmaMotor)),2); 
            likeV = lapseRV*lapsePDF + (1-lapseRV)*likeV;                                               %Correct the likelihoods for the lapse rate pdf
            LLvectorV(:) = log(likeV);                                                                  %Log transform and save

            % If requested as output: construct an approximate likelihood function following the alternative integration method (2. see below)
            if nargout > 3
                [likeFunV,likeFunV_grid] = BCIcompLikeFun(sV_hat,xV_prob,nGrid,RespRange,2,sigmaMotor);
                likeFunV = lapseRV*lapsePDF + (1-lapseRV)*likeFunV;                                     %Correct the likelihood function for the lapse rate
            end    

        % Alternative integration method for continuous responses    
        elseif integrateMethod == 2

            % Integrate out xV using a custom-made numerical intergration function that takes into account the width of each sV_hat bin in the xV direction    
            [likeFunV,likeFunV_grid,likeFunV_spacing] = BCIcompLikeFun(sV_hat,xV_prob,nGrid,RespRange,2,sigmaMotor);
            likeFunV = lapseRV*lapsePDF + (1-lapseRV)*likeFunV;                                         %Correct the likelihood function for the lapse rate
            LLvectorV(:) = log(lininterp1(likeFunV_grid,likeFunV,sV_resp,[],likeFunV_spacing));         %Interpolate likelihood function to find pdf value for each response 
        end
    end
    
    % Auditory responses (for unisensory visual stimuli) --> Set likelihood equal to lapse rate   
    if nTrialsA > 0 
        if discreteBool
            LLvectorA(:) = log(1/nRespLocs);
        else
            LLvectorA(:) = log(lapsePDF);
        end
    end   
    
    % Common source judgments (for unisensory visual stimuli)
    if nTrialsC > 0 
        likeFunCSJ = lapseRC*0.5 + (1-lapseRC)*[0 1];                                               %[Common = 0, Indep = 1];
        LLvectorC(:) = log(likeFunCSJ(CSJ_resp));                                                   %Save the log-likelihoods per response                             
    end                                                                                             %Note that CSJ_resp was already transformed to bin_idx numbers before
    
% Bisensory audiovisual stimulus (AV) 
else
    
    % Sample xA and xV on regular grids and compute the probability for each combination of xA,xV
    xA = linspace(sA-(gridCoverSD/2)*sigmaA,sA+(gridCoverSD/2)*sigmaA,nGrid)' + noiseShiftA;        %1st Dim (column vector)
    xV = linspace(sV-(gridCoverSD/2)*sigmaV,sV+(gridCoverSD/2)*sigmaV,nGrid) + noiseShiftV;         %2nd Dim (row vector)
    xA_prob = normpdf(xA,sA,sigmaA)*A_spacing; 
    xV_prob = normpdf(xV,sV,sigmaV)*V_spacing; 
    xAV_prob = bsxfun(@times,xA_prob,xV_prob);
    
    % Correct the likelihood distributions for their eccentricity-dependent bias
    xA = xA + deltaXA*(xA-muP);
    xV = xV + deltaXV*(xV-muP);
    
    % Bias by prior: compute maximum a posteriori (MAP) estimates given two independent causes (C=2)    
    sA_hat_indep = (xA/varA + muP/varP) * varA_pos;
    sV_hat_indep = (xV/varV + muP/varP) * varV_pos;
    
    % Perform multisensory integration and bias by prior: compute MAP estimates given one common cause (C=1)  
    sAV_hat_common = (bsxfun(@plus,xA/varA,xV/varV)+ muP/varP) * varAV_pos;

    % Compute log-likelihoods conditional on one common source or two independent sources: log() of P(xA,xV|C=1) or P(xA,xV|C=2)
    log_like_common = -(bsxfun(@minus,xV,xA).^2*varP + bsxfun(@plus,(xV-muP).^2*varA,(xA-muP).^2*varV)) / (2*varAV_PxGivenCommon) - log(2*pi*sqrt(varAV_PxGivenCommon));
    
    log_like_indepA = -((xA-muP).^2)/(2*varA_PxGivenIndep) - log(sqrt(2*pi*varA_PxGivenIndep));
    log_like_indepV = -((xV-muP).^2)/(2*varV_PxGivenIndep) - log(sqrt(2*pi*varV_PxGivenIndep));
    log_like_indep = bsxfun(@plus,log_like_indepA,log_like_indepV);
    
    % Compute posterior probability of common source: P(C=1|xA,xV)
    pC = exp( (log_like_common+log(Pcommon)) - logplusexp(log_like_common+log(Pcommon), log_like_indep+log(1-Pcommon)) );
    
    % Note how we avoided a multiplication/division of two small likelihoods in the above (c.f. that would have occured in the analogous code below and leads to numerical instability: e.g. pC = NaN or Inf)   
    % likelihood_common = exp(-(bsxfun(@minus,xV,xA).^2*varP + bsxfun(@plus,(xV-muP).^2*varA,(xA-muP).^2*varV)) / (2*varAV_PxGivenCommon)) / (2*pi*sqrt(varAV_PxGivenCommon));
    % likelihood_indep =  bsxfun(@times, exp(-((xA-muP).^2)/(2*varA_PxGivenIndep)) / sqrt(2*pi*varA_PxGivenIndep), exp(-((xV-muP).^2)/(2*varV_PxGivenIndep)) / sqrt(2*pi*varV_PxGivenIndep) );
    % pC = (likelihood_common*Pcommon)./(likelihood_common*Pcommon + likelihood_indep*(1-Pcommon));
    
    % Determine probabilities of common/segregated source judgements for non-Bayesian decision functions    
    if strcmpi(decisionFun,'ModelSelection')
        xAV_prob_common = (pC > 0.5).*xAV_prob;                             %Choose integration (or segregation) deterministically based on P(C=1|xA,xV) > 0.5
        xAV_prob_indep = (pC <= 0.5).*xAV_prob;
    elseif strcmpi(decisionFun,'ProbabilityMatching')
        xAV_prob_common = pC.*xAV_prob;                                     %Stochastically choose integration or segregation...
        xAV_prob_indep = (1-pC).*xAV_prob;                                  %where the probability for either choice is defined by P(C=1|xA,xV)
    end
    
    % Prepare some things once for discrete A and V responses...
    if discreteBool && ((nTrialsA > 0) || (nTrialsV > 0))
        % Prepare discrete response bins for integration in 3rd dimension
        Resp_edges_high = reshape(Resp_edges_high,[1 1 numel(Resp_edges_high)]);                                                    
        % Integrate out xV (2nd dim) and xA (1st dim) for common source responses only (this is not relevant for 'ModelAveraging') 
        if strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
            cdf_like_common = qtrapz(qtrapz(bsxfun(@times,xAV_prob_common,bsxfun(@ge,Resp_edges_high,sAV_hat_common)),2),1);   
        end                                                                                                                             
    end
        
    % Auditory responses
    if nTrialsA > 0

        % Compute final estimate for ModelAveraging decision function
        if strcmpi(decisionFun,'ModelAveraging')
            sA_hat = pC.*sAV_hat_common + bsxfun(@times,(1-pC),sA_hat_indep);                                   %Probability weighted average of integration and segregation estimates
        end 
        
        % Discrete responses
        if discreteBool
            
            % Apply trapezoidal integration over xA (1st dim) and xV (2nd dim) separately for each response bin (in the 3rd dimension). 
            if strcmpi(decisionFun,'ModelAveraging')
                cdf_likeA = qtrapz(qtrapz(bsxfun(@times,xAV_prob,bsxfun(@ge,Resp_edges_high,sA_hat)),2),1);     
            % Integrate out xA (1st dim) for segregated source judgments and sum likelihoods (per response option in 3rd dimension) over common and segregated source judgments (i.e. marginalize over the two options)  
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                cdf_likeA = cdf_like_common + qtrapz(bsxfun(@times,sum(xAV_prob_indep,2),bsxfun(@ge,Resp_edges_high,sA_hat_indep)),1);  
            end
            likeFunA = diff([0 reshape(cdf_likeA,[1 nRespLocs])]);                                              %Convert back to probabilities         
            likeFunA = lapseRA*(1/nRespLocs) + (1-lapseRA)*likeFunA;                                            %Correct the likelihood function for the lapse rate
            LLvectorA(:) = log(likeFunA(sA_resp));                                                              %Save the log-likelihoods per response                             
                                                                                                                %Note that sA_resp was already transformed to bin_idx numbers before
        % Default integration method for continuous responses     
        elseif integrateMethod == 1
            
            % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.  
            % For each "response bin" (with centres defined by sA_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sA_hat
            % We can then compute the overall likelihood for each response by integrating out the xVs (2nd dim) and xAs (1st dim) (i.e. the likelihood is a weighted average pdf of the motor noise Gaussian)   
            if strcmpi(decisionFun,'ModelAveraging')
                likeA = qtrapz(qtrapz(bsxfun(@times,xAV_prob,motorNoisePDFfun(reshape(sA_resp,[1 1 nTrialsA]), sA_hat, sigmaMotor)),2),1);
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                likeA_common = qtrapz(qtrapz(bsxfun(@times,xAV_prob_common,motorNoisePDFfun(reshape(sA_resp,[1 1 nTrialsA]), sAV_hat_common, sigmaMotor)),2),1);
                likeA_indep = qtrapz(bsxfun(@times,sum(xAV_prob_indep,2),motorNoisePDFfun(reshape(sA_resp,[1 1 nTrialsA]), sA_hat_indep, sigmaMotor)),1);
                likeA = likeA_common + likeA_indep;                                             %Sum likelihoods (per response) over common and segregated source judgments (i.e. marginalize over the two options)
            end
            LLvectorA(:) = log(lapseRA*lapsePDF + (1-lapseRA)*reshape(likeA,[nTrialsA 1]));     %Correct the likelihoods for the lapse rate pdf, log transform and save
            
            % If requested as output: construct an approximate likelihood function following the alternative integration method (2. see below)     
            if nargout > 3
                % Concatenate common and independent sources
                if strcmpi(decisionFun,'ModelAveraging')
                    xAV_prob_A = xAV_prob;
                elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                    sA_hat = cat(2,sAV_hat_common,sA_hat_indep);                                    
                    xAV_prob_A = cat(2,xAV_prob_common,sum(xAV_prob_indep,2));
                end
                % Integrate out xA and xV using a custom-made numerical intergration function that takes into account the width of each sA_hat bin in the xA direction (1st)   
                [likeFunA,likeFunA_grid] = BCIcompLikeFun(sA_hat,xAV_prob_A,nGrid,RespRange,1,sigmaMotor);
                likeFunA = lapseRA*lapsePDF + (1-lapseRA)*likeFunA;                                             %Correct the likelihood function for the lapse rate  
            end
            
        % Alternative integration method for continuous responses    
        elseif integrateMethod == 2
            
            % Concatenate common and independent sources
            if strcmpi(decisionFun,'ModelAveraging')
                xAV_prob_A = xAV_prob;
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                sA_hat = cat(2,sAV_hat_common,sA_hat_indep);                                    
                xAV_prob_A = cat(2,xAV_prob_common,sum(xAV_prob_indep,2));
            end
            % Integrate out xA and xV using a custom-made numerical intergration function that takes into account the width of each sA_hat bin in the xA direction (1st)   
            [likeFunA,likeFunA_grid,likeFunA_spacing] = BCIcompLikeFun(sA_hat,xAV_prob_A,nGrid,RespRange,1,sigmaMotor);
            likeFunA = lapseRA*lapsePDF + (1-lapseRA)*likeFunA;                                                 %Correct the likelihood function for the lapse rate
            LLvectorA(:) = log(lininterp1(likeFunA_grid,likeFunA,sA_resp,[],likeFunA_spacing));                 %Interpolate likelihood function to find pdf value for each response 
        end
    end
    
    % Visual responses
    if nTrialsV > 0

        % Compute final estimate for ModelAveraging decision function
        if strcmpi(decisionFun,'ModelAveraging')
            sV_hat = pC.*sAV_hat_common + bsxfun(@times,(1-pC),sV_hat_indep);                                   %Probability weighted average of integration and segregation estimates
        end
        
        % Discrete responses
        if discreteBool
            
            % Apply trapezoidal integration over xA (1st dim) and xV (2nd dim) separately for each response bin (in the 3rd dimension). 
            if strcmpi(decisionFun,'ModelAveraging')
                cdf_likeV = qtrapz(qtrapz(bsxfun(@times,xAV_prob,bsxfun(@ge,Resp_edges_high,sV_hat)),2),1);     
            % Integrate out xV (2nd dim) for segregated source judgments and sum likelihoods (per response option in 3rd dimension) over common and segregated source judgments (i.e. marginalize over the two options)  
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                cdf_likeV = cdf_like_common + qtrapz(bsxfun(@times,sum(xAV_prob_indep,1),bsxfun(@ge,Resp_edges_high,sV_hat_indep)),2);  
            end
            likeFunV = diff([0 reshape(cdf_likeV,[1 nRespLocs])]);                                              %Convert back to probabilities  
            likeFunV = lapseRV*(1/nRespLocs) + (1-lapseRV)*likeFunV;                                            %Correct the likelihood function for the lapse rate
            LLvectorV(:) = log(likeFunV(sV_resp));                                                              %Save the log-likelihoods per response                             
                                                                                                                %Note that sV_resp was already transformed to bin_idx numbers before   
        % Default integration method for continuous responses     
        elseif integrateMethod == 1
            
            % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.  
            % For each "response bin" (with centres defined by sV_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sV_hat
            % We can then compute the overall likelihood for each response by integrating out the xVs (2nd dim) and xAs (1st dim) (i.e. the likelihood is a weighted average pdf of the motor noise Gaussian)   
            if strcmpi(decisionFun,'ModelAveraging')
                likeV = qtrapz(qtrapz(bsxfun(@times,xAV_prob,motorNoisePDFfun(reshape(sV_resp,[1 1 nTrialsV]), sV_hat, sigmaMotor)),2),1);
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')                                      
                likeV_common = qtrapz(qtrapz(bsxfun(@times,xAV_prob_common,motorNoisePDFfun(reshape(sV_resp,[1 1 nTrialsV]), sAV_hat_common, sigmaMotor)),2),1);
                likeV_indep = qtrapz(bsxfun(@times,sum(xAV_prob_indep,1),motorNoisePDFfun(reshape(sV_resp,[1 1 nTrialsV]), sV_hat_indep, sigmaMotor)),2);
                likeV = likeV_common + likeV_indep;                                             %Sum likelihoods (per response) over common and segregated source judgments (i.e. marginalize over the two options)
            end
            LLvectorV(:) = log(lapseRV*lapsePDF + (1-lapseRV)*reshape(likeV,[nTrialsV 1]));     %Correct the likelihoods for the lapse rate pdf, log transform and save
            
            % If requested as output: construct an approximate likelihood function following the alternative integration method (2. see below)    
            if nargout > 3
                % Concatenate common and independent sources
                if strcmpi(decisionFun,'ModelAveraging')
                    xAV_prob_V = xAV_prob;
                elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                    sV_hat = cat(1,sAV_hat_common,sV_hat_indep);                                    
                    xAV_prob_V = cat(1,xAV_prob_common,sum(xAV_prob_indep,1));
                end
                % Integrate out xA and xV using a custom-made numerical intergration function that takes into account the width of each sV_hat bin in the xV direction (2nd)   
                [likeFunV,likeFunV_grid] = BCIcompLikeFun(sV_hat,xAV_prob_V,nGrid,RespRange,2,sigmaMotor);
                likeFunV = lapseRV*lapsePDF + (1-lapseRV)*likeFunV;                                             %Correct the likelihood function for the lapse rate
            end
                                                                                                                
        % Alternative integration method for continuous responses    
        elseif integrateMethod == 2
            
            % Concatenate common and independent sources
            if strcmpi(decisionFun,'ModelAveraging')
            	xAV_prob_V = xAV_prob;
            elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
                sV_hat = cat(1,sAV_hat_common,sV_hat_indep);                                    
                xAV_prob_V = cat(1,xAV_prob_common,sum(xAV_prob_indep,1));
            end
            % Integrate out xA and xV using a custom-made numerical intergration function that takes into account the width of each sV_hat bin in the xV direction (2nd)   
            [likeFunV,likeFunV_grid,likeFunV_spacing] = BCIcompLikeFun(sV_hat,xAV_prob_V,nGrid,RespRange,2,sigmaMotor);
            likeFunV = lapseRV*lapsePDF + (1-lapseRV)*likeFunV;                                             	%Correct the likelihood function for the lapse rate
            LLvectorV(:) = log(lininterp1(likeFunV_grid,likeFunV,sV_resp,[],likeFunV_spacing));                 %Interpolate likelihood function to find pdf value for each response 
        end
    end
    
    % Common source judgments
    if nTrialsC > 0
        if strcmpi(decisionFun,'ModelAveraging')
            likeFunCSJ(1) = qtrapz(qtrapz(xAV_prob.*(pC >= CSJthresh)));        %1 common source
            likeFunCSJ(2) = qtrapz(qtrapz(xAV_prob.*(pC < CSJthresh)));         %2 independent sources
        elseif strcmpi(decisionFun,'ModelSelection') || strcmpi(decisionFun,'ProbabilityMatching')
            likeFunCSJ(1) = qtrapz(qtrapz(xAV_prob_common));                    %1 common source
            likeFunCSJ(2) = qtrapz(qtrapz(xAV_prob_indep));                     %2 independent sources
        end
        likeFunCSJ = lapseRC*0.5 + (1-lapseRC)*likeFunCSJ;                      %Correct the likelihood function for the lapse rate
        LLvectorC(:) = log(likeFunCSJ(CSJ_resp));                               %Save the log-likelihoods per response                             
    end                                                                         %Note that CSJ_resp was already transformed to bin_idx numbers before
    
end %end if-statement: unisensory A, unisensory V, or bisensory AV location 

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function lpe = logplusexp(A, B)
% LOGPLUSEXP(A, B) computes log(exp(A) + exp(B)) robustly, actually log(bsxfun(@plus, exp(A), exp(B)))
%
%     lpe = logplusexp(A, B);
%
% Avoids overflow when both A and B are large.
%
% This function also avoids underflow when A or B is (close to) zero and the
% other is moderately negative. This case rarely comes up in my usage, and is
% neglected in the log(cum)sumexp routines.
%
% SEE ALSO: LOGSUMEXP LOGCUMSUMEXP

% Iain Murray, September 2010

mx = bsxfun(@max, A, B);
mn = bsxfun(@min, A, B);
%lcs = mx + log(1 + exp(mn-mx));
lpe = mx + log1p(exp(mn-mx));

% Fix up NaN and Inf special cases:
%inf_mask = bsxfun(@or, isinf(A), isinf(B));  %--> COMMENTED OUT FOR SPEED
%lpe(inf_mask) = mx(inf_mask);
%lpe(bsxfun(@or, isnan(A), isnan(B))) = NaN;

end %[EoF]

%-------------------------------------------------------------------------

function z = qtrapz(y,dim)
%QTRAPZ  Quick trapezoidal numerical integration.
%   Z = QTRAPZ(Y) computes an approximation of the integral of Y via
%   the trapezoidal method (with unit spacing).  To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, QTRAPZ(Y) is the integral of Y. For matrices, QTRAPZ(Y)
%   is a row vector with the integral over each column. For N-D
%   arrays, QTRAPZ(Y) works across the first non-singleton dimension.
%
%   Z = QTRAPZ(Y,DIM) integrates across dimension DIM of Y. The length of X 
%   must be the same as size(Y,DIM).
%
%   QTRAPZ is up to 3-4 times faster than TRAPZ for large arrays.
%
%   See also TRAPZ.

% Luigi Acerbi <luigi.acerbi@nyu.edu>
% Version 1.0. Release date: Jul/20/2015.

% By default integrate along the first non-singleton dimension
if nargin < 2; dim = find(size(y)~=1,1); end    

% Behaves as sum on empty array
if isempty(y); z = sum(y,dim); return; end

% Compute dimensions of input matrix    
if isvector(y); n = 1; else n = ndims(y); end

switch n
    case {1,2}      % 1-D or 2-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:) + y(end,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1) + y(:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end

    case 3      % 3-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:) + y(end,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:) + y(:,end,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1) + y(:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    case 4      % 4-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:,:) + y(end,:,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:,:) + y(:,end,:,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1,:) + y(:,:,end,:));
            case 4
                z = sum(y,4) - 0.5*(y(:,:,:,1) + y(:,:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    otherwise   % 5-D array or more
        for iDim = 1:n; index{iDim} = 1:size(y,iDim); end
        index1 = index;     index1{dim} = 1;
        indexend = index;   indexend{dim} = size(y,dim);
        try
            z = sum(y,dim) - 0.5*(y(index1{:}) + y(indexend{:}));
        catch
            error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');            
        end
end

end %[EoF]

%-------------------------------------------------------------------------

function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).
%   Y = BSXFUN_NORMPDF(X,MU,SIGMA) returns the pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMCDF, NORMPDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
    elseif isscalar(sigma)
        y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
    else
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

end %[EoF]

%-------------------------------------------------------------------------

function y = bsxfun_triangpdf(x,mu,sigma)
%BSXFUN_TRIANGPDF Vectorized probability density function (pdf) of
%symmetrical triangular distribution as an approximation to the normal pdf. 
%   Y = BSXFUN_TRIANGPDF(X,MU,SIGMA) returns the pdf of the symmetrical
%   triangular distribution with mean MU. Its width is defined relative to 
%   width parameter SIGMA such that the interval MU +/- SIGMA contains ~68%
%   of its mass (similar to a normal distribution with same mu and sigma).
%   Dimensions of X, MU, and SIGMA must either match, or be equal to one. 
%   Computation of the pdf is performed with singleton expansion enabled 
%   via BSXFUN. The size of Y is the size of the input arguments (expanded
%   to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   Adapted from BSXFUN_NORMPDF by Luigi Acerbi.
%
%   This function is about 3 times faster than BSXFUN_NORMPDF

if nargin<3
    error('bmp:bsxfun_triangpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

widthFactor = 2.29; %--> Create a symmetrical triangular pdf where 68% of its mass falls within the interval that would also contain 68% of the normal distribution's mass (i.e. interval is mu +/- SD)
                         %This also means that 97.8% of the gaussian mass is captured within the interval of the triangular pdf (outside that interval pdf = 0 for triangular, and pdf = small for gaussian)   

%y = ((widthFactor*sigma)-abs(mu-x))/((widthFactor*sigma)^2)
try
    if isscalar(mu)
        y = bsxfun(@rdivide,bsxfun(@minus,(widthFactor*sigma),abs(mu-x)),(widthFactor*sigma).^2);        
    elseif isscalar(sigma)
        y = ((widthFactor*sigma)-abs(bsxfun(@minus,mu,x)))/((widthFactor*sigma)^2);
    else
        y = bsxfun(@rdivide,bsxfun(@minus,(widthFactor*sigma),abs(bsxfun(@minus,mu,x))),(widthFactor*sigma).^2);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

y = max(y,0);       %This is significantly (!) faster than y(y<0) = 0;

end %[EoF]

%-------------------------------------------------------------------------

function Vout = lininterp1(X,V,Xq,extrap,deltaX)
%Simple 1D linear interpolation on a regular ascending grid.
%Inputs X, V and Xq must be vectors. 
%Output Vout is a column vector.

%This function is a simplified version of the lininterp1 function that was
%originally written by Luigi Acerbi: https://github.com/lacerbi/lautils-mat
%While Luigi's version supports matrices as inputs, this simplified
%version is about 4x faster, and nearly 10x faster than Matlab's interp1.

%Determine length of the given grid 
%Note that X is allowed to be given as X = [min max];
nGrid = length(V);

%Ensure Xq is not empty and a column vector
if isempty(Xq)
    Vout = [];
    return
else
    Xq = Xq(:);
end

%Set default extrapolation: NaN
if nargin < 4
    extrap = NaN;
end

%Determine deltaX from X if not already supplied
if nargin < 5
    deltaX = (X(end)-X(1))/(nGrid-1);
end

%Find Xq out of bounds
flag1 = Xq < X(1);
flag2 = Xq > X(end);
flag3 = isnan(Xq);

%Transform Xq as if X was given as 1:length(X)
Xq = (Xq-X(1))/deltaX + 1;                                                                          
                                
%Indices of nearest lower neighbour X
Xqi = floor(Xq);
Xqi(flag1 | flag2 | flag3) = 1;     %temp --> will be overwritten below

%Distance from Xqi
delta = Xq - Xqi; 

%Expand to avoid Xqi==numel(V) errors when asking for Xqi+1 index below
V = [V(:); 0];

%Distance-weighted average (i.e. linear interpolation)   
Vout = (1-delta).*V(Xqi) + delta.*V(Xqi+1);

%Set values outside of X (extrapolation)
if isempty(extrap)                                              
    %Xq outside X get V values of nearest edge.
    Vout(flag1) = V(1);        
    Vout(flag2) = V(nGrid);
else
    %Xq outside X get specified value (default = NaN)
    Vout(flag1 | flag2) = extrap;
end
Vout(flag3) = NaN; 

end %[EoF]
