function mcmcStructure = eissample_wrapper(logPfuns,x0,N,K,widths,LB,UB,options,nAttempts,mcmcStructureLastTime,varargin)
%Small wrapper around EISSAMPLE in which we attempt to reach convergence
%by adaptively increasing thinning. We DO NOT make use of parallelization.

%You can also continue a previous sampling action for this same dataset:
%simply increase the value of "nAttempts" before calling this function 
%again and add the previous results in "mcmcStructureLastTime" (empty 
%otherwise; don't change any of the other settings or it may not work).

%Shuffle the random number generator (in case this function is called in parallel / multiple times - just ensure some time difference between each succesive/parallel call)      
rng('shuffle');                                                                                                                 

%Initialize output
mcmcStructure = [];

%Unpack nAttempts
minRounds = nAttempts(1);
maxRounds = nAttempts(2);

%Ensure that EISSAMPLE does not perform the diagnostics (we do it ourselves in this function)   
options.Diagnostics = false;

%Throw a warning if thinning is requested in EISSAMPLE and adaptive thinning is requested by calling this function   
if isfield(options,'Thin')
    if (options.Thin > 1) && (maxRounds > 1)
        warning('Thinning is requested to EISSAMPLE but we apply additional thinning if MCMC does not converge in the first round. Consider setting "options.Thin = 0" and let this function handle thinning instead.');
    end
end

%Set some defaults, ensure these fields exist
if ~isfield(options,'VPbasedStepOut')
    options.VPbasedStepOut = false;
end
if ~isfield(options,'Noise')
    options.Noise = false;  %This one is used for calling "samplediagnose_DM"
end
if isfield(options,'Display')
    switch options.Display
        case {'iter','iter-detailed'}
            trace = 3;
        case {'notify','notify-detailed'}
            trace = 2;
        case {'final','final-detailed'}
            trace = 1;
        case {'none', 'off'}
            trace = 0;
        otherwise
            trace = 1;
    end
else
    trace = 2;
end
if ~iscell(logPfuns)
    logPfuns = {logPfuns};  %We assume this further below
    options.nocell = 1;     %This option is used to unpack fvals in eissample_lite
end

%Correct the number of samples such that it is a multiple of nWalkers
samplesPerWalker = ceil(N/K);
mccountRounded = K*samplesPerWalker;  

%Try to obtain the posterior using MCMC, repeat if necessary
nNewSamplesPerRound = [1 2.^(0:(maxRounds-2))];                             %i.e. the number of times that we call the eissample function per round is: [1 1 2 4 8 etc.]   
thinningPerRound = 2.^(0:(maxRounds-1));                                    %i.e. thinning per round is: [1 2 4 8 16 etc.]

%Should we continue with a previous set of MCMC samples and sample more rounds with more thinning?   
if ~isempty(mcmcStructureLastTime)
    
    %Unpack old results
    mcmcParams = mcmcStructureLastTime.mcmcParams;
    mcmcLP = mcmcStructureLastTime.prob.mcmcLogPosterior;
    mcmcLPrior = mcmcStructureLastTime.prob.mcmcLogPrior;
    mcmcLL_responses = mcmcStructureLastTime.prob.mcmcLL_responses;
    
    %Ensure that we don't repeat the Burn-in
    BurnInCopy = mcmcStructureLastTime.output.burn;
    options.Burnin = K;
    
    %Overwrite some other settings
    if ~options.VPbasedStepOut
        widths = mcmcStructureLastTime.output.covmat;                       %The covmat may have been updated adaptively during a previous burn-in.
    end
    output = mcmcStructureLastTime.output;
    
    %Initialize some counters
    MCMCroundCounter = mcmcStructureLastTime.output.nRounds + 1;
    funccount = mcmcStructureLastTime.output.funccount;
    nslicecollapsed = mcmcStructureLastTime.output.nslicecollapsed;
    nTotSamples = mcmcStructureLastTime.output.totN;    
    SigmaFactor = mcmcStructureLastTime.output.SigmaFactor; 
    SigmaFactorCount = mcmcStructureLastTime.output.SigmaFactorCount; 
else
    %Initialize some counters
    MCMCroundCounter = 1;
    funccount = zeros(1,numel(logPfuns));                                   %Note that we assume that the logPfuns are put into cell-arrays
    nslicecollapsed = 0;
    nTotSamples = 0;
    SigmaFactor = 0;
    SigmaFactorCount = 0;
end
MCMCconvergedBool = 0;  %Even if it has already converged last time, the fact that this function was called again probably means that more samples should be acquired...

%Start main while loop: Collect samples using MCMC!
while ~MCMCconvergedBool && (MCMCroundCounter <= maxRounds)
    
    %Start a new timer on every round
    cStartRound = clock;

    %Initialize
    if MCMCroundCounter == 1
        mcmcParams_Cell = cell(1,K);
        mcmcLP_Cell = cell(1,K);
        mcmcLPrior_Cell = cell(1,K);
        mcmcLL_responses_Cell = cell(1,K);
    else
        %Unmerge the chains from the previous round to form cell arrays of size [1 x K]
        mcmcParams_Cell = mergeChains(mcmcParams,N,K,1);
        mcmcLP_Cell = mergeChains(mcmcLP,N,K,1);
        mcmcLPrior_Cell = mergeChains(mcmcLPrior,N,K,1);
        mcmcLL_responses_Cell = mergeChains(mcmcLL_responses,N,K,1);
        clear mcmcParams mcmcLP mcmcLPrior mcmcLL_responses
        
        %Create new initializations for the walkers on this new round
        x0 = nan(K,size(mcmcParams_Cell{1},2));
        for k=1:K
            x0(k,:) = mcmcParams_Cell{k}(end,:);
        end
        
        %Thin each chain/walker by half (the thinned/deleted samples will be replaced with newly collected samples this round)    
        for k = 1:K
            mcmcParams_Cell{k} = mcmcParams_Cell{k}(1:2:end,:); 
            mcmcLP_Cell{k} = mcmcLP_Cell{k}(1:2:end,:); 
            mcmcLPrior_Cell{k} = mcmcLPrior_Cell{k}(1:2:end,:); 
            mcmcLL_responses_Cell{k} = mcmcLL_responses_Cell{k}(1:2:end,:); 
        end
    end 
    
    %Sample new MCMC samples and take care of thinning
    for i_thin=1:nNewSamplesPerRound(MCMCroundCounter)
        
        %Call EISSAMPLE (by Luigi Acerbi, https://github.com/lacerbi/eissample/). The function was adjusted by David Meijer
        [mcmcParams_Tmp,mcmcLP_Tmp,~,output_Tmp,mcmcLL_responses_temp] = eissample_lite_DM2(logPfuns,x0,mccountRounded,K,widths,LB,UB,options,varargin{:});
        mcmcLPrior_Tmp = mcmcLL_responses_temp{1};            
        mcmcLL_responses_Tmp = mcmcLL_responses_temp{2};             
            
        %Create an output structure
        if (MCMCroundCounter == 1)
            output = output_Tmp;
            BurnInCopy = output_Tmp.burn;                                   %Reset the burnin options setting (avoid using a Burnin period on consecutive calls of EISSAMPLE)
            options.Burnin = K;                                             %The first sample of every walker is burned - these were the last samples of the previous round..
        end                                                                 
        
        %Update the function counts  
        funccount = funccount + output_Tmp.funccount;                      
        nslicecollapsed = nslicecollapsed + output_Tmp.nslicecollapsed;
        nTotSamples = nTotSamples + output_Tmp.totN;
        CountTemp = sum(~isnan(output_Tmp.SigmaFactor));
        SigmaFactor = (SigmaFactorCount*SigmaFactor + CountTemp*nanmean(output_Tmp.SigmaFactor)) / (SigmaFactorCount + CountTemp);
        SigmaFactorCount = SigmaFactorCount + CountTemp;
        
        %Split the results into its individual chains/walkers (i.e. unmerge chains)    
        mcmcParams_Tmp_Cell = mergeChains(mcmcParams_Tmp,size(mcmcParams_Tmp,1),K,1);
        mcmcLP_Tmp_Cell = mergeChains(mcmcLP_Tmp,size(mcmcLP_Tmp,1),K,1);
        mcmcLPrior_Tmp_Cell = mergeChains(mcmcLPrior_Tmp,size(mcmcLPrior_Tmp,1),K,1);
        mcmcLL_responses_Tmp_Cell = mergeChains(mcmcLL_responses_Tmp,size(mcmcLL_responses_Tmp,1),K,1);
        clear mcmcParams_Tmp mcmcLP_Tmp mcmcLPrior_Tmp mcmcLL_responses_Tmp
        
        %Thin and save (simply concatenate)
        ThinBy = thinningPerRound(MCMCroundCounter);
        for k = 1:K 
            mcmcParams_Cell{k} = [mcmcParams_Cell{k}; mcmcParams_Tmp_Cell{k}(1:ThinBy:end,:)]; 
            mcmcLP_Cell{k} = [mcmcLP_Cell{k}; mcmcLP_Tmp_Cell{k}(1:ThinBy:end,:)]; 
            mcmcLPrior_Cell{k} = [mcmcLPrior_Cell{k}; mcmcLPrior_Tmp_Cell{k}(1:ThinBy:end,:)];
            mcmcLL_responses_Cell{k} = [mcmcLL_responses_Cell{k}; mcmcLL_responses_Tmp_Cell{k}(1:ThinBy:end,:)]; 
        end
        clear mcmcParams_Tmp_Cell mcmcLP_Tmp_Cell mcmcLPrior_Tmp_Cell mcmcLL_responses_Tmp_Cell
    
    end %End of 1 round of MCMC
    
    %Merge the chains to form simple matrices
    mcmcParams = mergeChains(mcmcParams_Cell,N,K);
    mcmcLP = mergeChains(mcmcLP_Cell,N,K);
    mcmcLPrior = mergeChains(mcmcLPrior_Cell,N,K);
    mcmcLL_responses = mergeChains(mcmcLL_responses_Cell,N,K);
    clear mcmcParams_Cell mcmcLP_Cell mcmcLPrior_Cell mcmcLL_responses_Cell
    
    %Check convergence: i.e. compute 'potential scale reduction factor' (R) and effective sample size (Neff) --> see Gelman et al., 2013 (book)
    [exitflag,R,Neff,tau] = samplediagnose_DM(mcmcParams,K,funccount,nslicecollapsed,trace,options.Noise);
    if trace > 0; disp(' '); end
    
    %Report timing of this round
    fprintf('Finished round %d, elapsed time (days hours:minutes:seconds) %s \n',MCMCroundCounter,datestr(etime(clock,cStartRound)/86400,'dd HH:MM:SS'));

    %Converged?
    if (exitflag == 1) && (MCMCroundCounter >= minRounds)
        MCMCconvergedBool = 1;                  %One way to leave the while loop
    %New round?
    else
        if MCMCroundCounter < maxRounds
            disp(['MCMC did not converge in round ' num2str(MCMCroundCounter) '. Starting a new round...']);
        else
            disp('MCMC did not converge in either of the rounds! Check the data and try again later using different settings (e.g. increase "nAttempts.MCMC").');
        end
        MCMCroundCounter = MCMCroundCounter+1;  %Potential other way to leave the while loop
    end 
end %End of while loop: converged? or maxRuns reached?

%Subtract one counter if we exited the while loop because MCMCroundCounter > maxRounds (new round was not actually executed)   
if ~MCMCconvergedBool
    MCMCroundCounter = MCMCroundCounter-1;
end

%Check that eissample_lite did run
if ~exist('exitflag','var')
    warning('Eissample did not run! Presumably because "MCMCroundCounter > maxRounds". Increase maxRounds to continue running an old MCMC sequence.');
    [exitflag,R,Neff,tau] = samplediagnose_DM(mcmcParams,K,funccount,nslicecollapsed,trace,options.Noise);  %Ensure that all the necessary output variables exist
end

%Update the output structure
output.nRounds = MCMCroundCounter;
output.burn = BurnInCopy;
output.thin = thinningPerRound(MCMCroundCounter);
output.N = N;
output.totN = nTotSamples;
output.funccount = funccount;
output.nslicecollapsed = nslicecollapsed;
output.SigmaFactor = SigmaFactor;
output.SigmaFactorCount = SigmaFactorCount;
output.R = R;
output.Neff = Neff;
output.tau = tau;

%Collect all variables in struct 'mcmcStructure' for the output of this function
mcmcStructure.mcmcParams = mcmcParams;

mcmcStructure.prob.mcmcLogPrior = mcmcLPrior;
mcmcStructure.prob.mcmcLogLikelihood = mcmcLP-mcmcLPrior;                   %Should be equal to "sum(mcmcLL_responses,2)"
mcmcStructure.prob.mcmcLogPosterior = mcmcLP;
mcmcStructure.prob.mcmcLL_responses = mcmcLL_responses;

mcmcStructure.converge.exitflag = exitflag;
mcmcStructure.converge.maxR = max(R);
mcmcStructure.converge.minNeff = min(Neff);

mcmcStructure.output = output;

% %Plot the results
% bounds = quantile(mcmcStructure.mcmcParams,[0.05 0.95],1);                %posterior quantiles (0.5 = median). This returns a matrix with size [2 x nParams]
% cornerplot_DM(mcmcStructure.mcmcParams,[],[],bounds,[],mcmcStructure.prob.mcmcLogPosterior);

end %[EOF]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function mergedData_Matrix = mergeChains(DataPerWalker_Cell,mccount,nWalkers,unMergeBool)
%Merge various walkers of MCMC sampler. Ensure that the number of samples is exactly mccount. 
%This function can also be used to "unmerge" the samples per Walker (set: unMergeBool = 1)

%Set default direction: "merge"
if nargin < 4
    unMergeBool = 0;
end

%Prepare to merge
minSamplesPerWalker = floor(mccount/nWalkers);
leftOver = mccount-(minSamplesPerWalker*nWalkers);
SamplesEachWalker = minSamplesPerWalker*ones(1,nWalkers);
SamplesEachWalker(1:leftOver) = SamplesEachWalker(1:leftOver)+1;            %The first few walkers get to keep 1 more sample than the last few walkers

%Merge
if unMergeBool == 0
    mergedData_Matrix = nan(mccount,size(DataPerWalker_Cell{1},2));
    for k=1:nWalkers
        idx = k:nWalkers:mccount;
        mergedData_Matrix(idx,:) = DataPerWalker_Cell{k}(1:SamplesEachWalker(k),:);
    end
%Unmerge    
elseif unMergeBool == 1
    mergedData_Matrix = cell(1,nWalkers);
    for k=1:nWalkers
        mergedData_Matrix{k} = DataPerWalker_Cell(k:nWalkers:mccount,:);
    end
end

end %[EoF]
