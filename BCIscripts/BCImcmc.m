function dataFit_MCMC = BCImcmc(BCIfitResults)
%Perform Markov Chain Monte Carlo to obatin an estimate of the
%posterior probability distribution over parameters.

%Shuffle the random number generator (ensure random initializations on multiple parallel calls of this function)     
rng('shuffle');

%Start a timer
disp('Starting MCMC fit ...');
cStart = clock;  

%Add relevant 'MCMC' subfolders to the path
me = mfilename;                                                 %what is my filename
pathstr = fileparts(fileparts(which(me)));                      %get location of 'BCImodel' folder
addpath(genpath([pathstr filesep 'eissample_lite_DM2']));         
addpath(genpath([pathstr filesep 'vbmc']));                     %Necessary for "VPbasedStepOut"
addpath(genpath([pathstr filesep 'psisloo_DM']));         
addpath(genpath([pathstr filesep 'gofit']));                    %Goodness-of-fit based on PSISLOO

%Get the data and some settings
[P,trueLocsAV,responsesAVC,i_Conditions] = unpackBCIfitResults(BCIfitResults);
LB = P.Bounds.conv.LB;
UB = P.Bounds.conv.UB;

%Define the functions that compute the log-likelihood and log-prior
LLfun = @(params,refLL) BCIcompLL(params,P,trueLocsAV,responsesAVC,i_Conditions,refLL);         %Returns a log-likelihood vector (one LL value for each response)
LPriorfun = @(params) BCIcomputeLogPrior(params,P);                                             %Returns the log-prior probability
logPfuns = {LPriorfun,LLfun};                                                                   %Note, the first function should be the prior (that one is evaluated first by eissample_lite.m)

%Some settings for the MCMC algorithm
if isfield(P,'MCMCcount')                                           %This is the number of MCMC samples returned after burn-in and after thinning
    mccount = P.MCMCcount;
else
    mccount = 10000;
end
nWalkers = 2*(P.nParams2Fit+1);                                     %This is the default for Eissample
options.Burnin = 0;                                                 %By default we don't use a Burnin/WarmUp period. However, if VBMC was not performed first, then we run an initial MCMC round with BurnIn (see below)
options.Thin = 1;                                                   %We don't use thinning. However, we'll use our own thinning implementation if MCMC did not converge the first time (see wrapper function) 
options.Display = 'off';                                            %No progress in Command Window please
options.Diagnostics = false;                                        %Do not let Eissample perform convergence checks. Instead we compute them ourselves (see wrapper function). We do this because of our thinning procedure..
options.Noise = true;                                               %We falsely claim that the LL function is stochastic. This doesn't actually change anything in the algortihm; it merely suppresses some warnings.
options.UseRefLogP = [false true];                                  %The log-prior function does not accept a logP reference value as second input parameter, but the LL function does.

%Get an estimate of the widths of the parameter distributions (this is preferable for the MCMC algorithm)
if isfield(BCIfitResults,'VBMC')
    options.FixedSliceLengths = true;                               %Use the covariance matrix (as defined by the VP) to fix the slice lengths (rather than using the default random slice lengths)
%     if BCIfitResults.VBMC.validation.exitflag == 1
%         options.VPbasedStepOut = true;                            %If the VBMC solution has been validated, use that VP to guide in setting the lengths of the slices and so reduce the number of inefficient function calls
%         widths = BCIfitResults.VBMC.vp;
%         x0 = vbmc_rnd(BCIfitResults.VBMC.vp,nWalkers,1,1);        %Sample initializations from the validated variational posterior (VBMC)
%     else
        options.VPbasedStepOut = false;                             %Alternatively, if the variational posterior was not validated (usually) ...
        [widths,x0] = compCovVBMC(BCIfitResults.VBMC,nWalkers,inf); %it often happens that one parameter does not contribute much to the LL, and therefore its value/distribution fluctuates across several VBMC runs.   
%     end                                                           %So, we compute covariance matrix and select initializations from multiple converged VPs (weighted by their probabilities / elbo) - see helper function below 
else
    %No VBMC has been performed
    options.Burnin = mccount;                                       %Run a warmp-up round to estimate correct widths
    options.FixedSliceLengths = true;                               %After the warm-up we can trust the covariance matrix to fix the slice lengths
    
    %Set some initial widths and samples X0
    fittedParamsBads = BCIconvertParams2Fit(BCIfitResults.BADS.fittedParams,P,'real2conv');
    [largerPLB,largerPUB] = shrinkPlausibleBounds(P,fittedParamsBads,10);           %This one is used for setting the widths of the parameters (but they will be reset adaptively during BurnIn)
    widths = largerPUB-largerPLB;
    
    [smallerPLB,smallerPUB] = shrinkPlausibleBounds(P,fittedParamsBads,100);        %This one is used to sample MCMC starting positions (initializations)
    x0 = sampleXnormal(fittedParamsBads,nWalkers,smallerPLB,smallerPUB);
end

%Should we continue with a previous set of MCMC samples and sample more rounds with more thinning?   
if isfield(BCIfitResults,'MCMC')
    mcmcStructureLastTime = BCIfitResults.MCMC;
    mcmcStructureLastTime.mcmcParams = BCIconvertParams2Fit(BCIfitResults.MCMC.mcmcParams,P,'real2conv');   %Convert variables to 'converted space'
else
    mcmcStructureLastTime = [];
end

%Call EISSAMPLE wrapper function --> perform MCMC (may take quite a while...)
if P.parallel.MCMC
    dataFit_MCMC = eissample_wrapper_par(logPfuns,x0,mccount,nWalkers,widths,LB,UB,options,P.nAttempts.MCMC,mcmcStructureLastTime); %Invalid method but faster convergence (assumes a valid VBMC estimate!)
    dataFit_MCMC.output.invalidParallelizationSpeedUp = 1;
else
    dataFit_MCMC = eissample_wrapper(logPfuns,x0,mccount,nWalkers,widths,LB,UB,options,P.nAttempts.MCMC,mcmcStructureLastTime);     %Default
end                                                                                                                                 

%Recompute Goodness-of-Fit using PSISLOO cross-validated log-likelihoods (use original settings and conditions structure - not the collapsed one, see "unpackBCIfitResults")  
%We use LLs (not posterior probabilities) as input to compute predictive densities because "we are interested here in summarizing the fit of model to data, and for this purpose 
%the prior is relevant in estimating the parameters [distribution] (i.e. estimating the posterior using MCMC) but not in assessing a model’s accuracy." (Gelman, Hwang, Vehtari, 2013)     
dataFit_MCMC.goodnessOfFit = BCIcomputeGoF(BCIfitResults.settings,trueLocsAV,responsesAVC,BCIfitResults.data.i_Conditions,'MCMC', ...
                                           BCIfitResults.BADS.goodnessOfFit,dataFit_MCMC.prob.mcmcLL_responses,dataFit_MCMC.converge.minNeff);

%Convert variables back to real space
dataFit_MCMC.mcmcParams = BCIconvertParams2Fit(dataFit_MCMC.mcmcParams,P,'conv2real');

%Report computation time in command window
fprintf('Finished MCMC fit, elapsed time (days hours:minutes:seconds) %s \n',datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end %[EOF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function [P,trueLocsAV,responsesAVC,i_Conditions] = unpackBCIfitResults(BCIfitResults)

%Get some settings
P = BCIfitResults.settings;

%Get the data
trueLocsAV = BCIfitResults.data.trueLocsAV;
responsesAVC = binResponses(P,BCIfitResults.data.responsesAVC);
i_Conditions = BCIfitResults.data.i_Conditions;

%Check whether certain conditions code for exactly the same parameters. 
%If so, merge the conditions. This is a significant speed-up!
ParamsPerCond = P.ParamsPerCond;
nCond = length(ParamsPerCond);
Conds2stay = 1:nCond;
for i=nCond:-1:2
    for j=(i-1):-1:1
        if isequal(sort(ParamsPerCond{i}),sort(ParamsPerCond{j}))
            i_Conditions(:,j) = i_Conditions(:,j) | i_Conditions(:,i);      %Add trials of deleted condition to other condition that fits exactly the same parameters (i.e. effectively that is no other condition)
            Conds2stay(i) = [];                                             %Delete duplicate condition
            break;                                                          %Break from inner for-loop
        end
    end
end
i_Conditions = i_Conditions(:,Conds2stay);
P.ParamsPerCond = ParamsPerCond(Conds2stay);
P.nConditions = length(Conds2stay);

end %[EoF]

%-------------------------------------------------------------------------

function binIdxResponsesAVC = binResponses(P,responsesAVC)

%Discrete responses
if all(~isnan(P.RespLocs))
    if numel(P.RespLocs) == 1
        Resp_edges = [-inf inf]; 
    else
        Resp_edges = [-inf mean([P.RespLocs(1:(end-1)); P.RespLocs(2:end)],1) inf];  
    end
    binIdxResponsesAVC = zeros(size(responsesAVC));
    [~,~,binIdxResponsesAVC(:,1)] = histcounts(responsesAVC(:,1),Resp_edges);
    [~,~,binIdxResponsesAVC(:,2)] = histcounts(responsesAVC(:,2),Resp_edges);
    binIdxResponsesAVC((responsesAVC(:,3) == 1),3) = 1;
    binIdxResponsesAVC((responsesAVC(:,3) ~= 1) & ~isnan(responsesAVC(:,3)),3) = 2;
    
    %histcounts returns zero for NaNs - Transform those back to NaNs
    binIdxResponsesAVC(binIdxResponsesAVC == 0) = NaN;
    
%Continuous responses    
else 
    binIdxResponsesAVC = responsesAVC;
    
    binIdxResponsesAVC((responsesAVC(:,3) == 1),3) = 1;                                 %A "1" means common-source
    binIdxResponsesAVC((responsesAVC(:,3) ~= 1) & ~isnan(responsesAVC(:,3)),3) = 2;     %Any other (e.g. "0" or "2") but not NaN means independent sources!
end

end %[EoF]

%-------------------------------------------------------------------------

function [newPLB,newPUB] = shrinkPlausibleBounds(P,meanX,factorSmaller)
%Decrease the space within the plausible bounds by some factorSmaller 
%Use meanX as centre point for the new plausible bounds. 
%Also ensure that meanX is within the hard bounds.

%Get the current bounds
LB = P.Bounds.conv.LB;
UB = P.Bounds.conv.UB;
PLB = P.Bounds.conv.PLB;
PUB = P.Bounds.conv.PUB;

%Size of current/old plausible bounds
rangeOld = PUB-PLB;
rangeNew = rangeOld/factorSmaller;
rangeEitherSide = rangeNew/2;

%Check that meanX is within the hard bounds (i.e. not ON the hard bounds),
%if on/outside, move meanX to halfway the hard and current/old plausible bounds    
meanX(meanX <= LB) = mean([LB(meanX <= LB); PLB(meanX <= LB)],1);           
meanX(meanX >= UB) = mean([LB(meanX >= UB); PLB(meanX >= UB)],1);           %Note that we assumed that [meanX,PLB,PUB,LB,UB] are all row vectors     

%Check that meanX is within the current/old plausible bounds, 
%if outside, move the current/old plausible bounds to halfway meanX and the hard bounds
PLB(meanX < PLB) = mean([meanX(meanX < PLB); LB(meanX < PLB)],1);           
PUB(meanX > PUB) = mean([meanX(meanX > PUB); UB(meanX > PUB)],1);           %Note that we assumed that [meanX,PLB,PUB,LB,UB] are all row vectors     

%Suggest new plausible bounds (we can only make the plausible space smaller, not larger)
newPLB = max([PLB; meanX-rangeEitherSide],[],1);
newPUB = min([PUB; meanX+rangeEitherSide],[],1);                            %Note that we assumed that [meanX,PLB,PUB] are all row vectors 

end %[EOF]

%-------------------------------------------------------------------------

function X = sampleXnormal(meanX,N,PLB,PUB)
%Sample from multivariate normal distribution with SDs set as: [PUB-PLB]/6. 
%Ensure that all samples fall within PUB and PLB

%Check that meanX is within the bounds
if any(meanX < PLB) || any(meanX > PUB)
    error('Error when creating intializations for MCMC: meanX falls outside of bounds [PLB PUB]');
end

%Create covariance matrix
SDs = (PUB-PLB)/6;
covmat = diag(SDs.^2);

%Sample until all X fall within bounds
X = nan(N,size(meanX,2));
invalid = true(N,1);
maxAttempts = 1e5;
attemptCounter = 0;
while any(invalid) && (attemptCounter < maxAttempts)
    attemptCounter = attemptCounter+1;
    idx_invalid = find(invalid);
    X(idx_invalid,:) = mvnrnd(meanX,covmat,numel(idx_invalid));  
    invalid = any(X <= PLB,2) | any(X >= PUB,2);
end

%Throw error if not all valid
if any(invalid)
    error('Error when creating intializations for MCMC: unable to find X within bounds [PLB PUB]');
end

end %[EoF] 

%-------------------------------------------------------------------------

function [covmat,x0] = compCovVBMC(vbmcStructure,nWalkers,useNrOfVPs)

%Set number of samples
nS = 1e5;

%Use only one VP (the best)
if useNrOfVPs == 1
    Xs = vbmc_rnd(vbmcStructure.vp,nS,1,1);                                 %Generate nS samples from the variational posterior
else
    %Check the number of available VPs
    nVPs = vbmcStructure.wrapper.nFinished;
    useNrOfVPs = min(nVPs,useNrOfVPs);
    
    %Check which VPs converged
    i_converged = vbmcStructure.wrapper.exitflag == 1;
    if sum(i_converged) == 0 
        i_converged = true(size(i_converged));                              %If none converged, than use all VPs (not recommended!). 
    end

    %Find the "useNrOfVPs" best VPs in terms of their lower confidence bound on ELBO
    elcbo = vbmcStructure.validation.elcbo;
    elcbo(~i_converged) = -Inf;
    [elcbo_sorted,idx_sort] = sort(elcbo,'descend');
    idx2use = idx_sort(1:useNrOfVPs);
    elcbo2use = elcbo_sorted(1:useNrOfVPs);
    
    %Sample a probability weighted number of samples from all chosen VPs
    weights = exp(elcbo2use-max(elcbo2use));                                %Convert ELBO to probabilities (p=1 for the highest elcbo)
    Xs = cell(1,useNrOfVPs);                    
    for i=1:useNrOfVPs
        vp_temp = vbmcStructure.wrapper.vp{idx2use(i),1};
        nS_temp = round(weights(i)*nS);
        if nS_temp > 0
            Xs{1,i} = vbmc_rnd(vp_temp,nS_temp,1,1);                        %Generate nS_temp random samples from the variational posterior   
        end
    end

    %Concatenate samples across VPs
    Xs = cat(1,Xs{:});
end

%Compute the covariance matrix
covmat = cov(Xs);                                     

%Randomly select initializations for the walkers
totNs = size(Xs,1);
idx_chosen = randperm(totNs,nWalkers);
x0 = Xs(idx_chosen,:);

end %[EoF]
