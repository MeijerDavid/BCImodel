function mcmcStructure = eissample_wrapper_par(logPfuns,x0,N,K,widths,LB,UB,options,nAttempts,mcmcStructureLastTime,varargin)
%Small wrapper around EISSAMPLE in which we attempt to reach convergence
%by adaptively increasing thinning. We make use of parallel computing by
%initializing walkers separately for each CPU worker. This method breaks
%the Markov Chains and is therefore invalid. However, if the VP result can
%be trusted (i.e. validated), it shouldn't affect the results, and it has 
%the advantage of faster convergence of the chains!

%Burnin is first distributed across all parallel workers.
%You can also continue a previous sampling action for this same dataset:
%simply increase the value of "maxRounds" before calling this function 
%again and add the previous results in "mcmcStructureLastTime" (empty 
%otherwise; don't change any of the other settings or it may not work).

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
if ~isfield(options,'Burnin')
    options.Burnin = 0;
end
if ~isfield(options,'Noise')
    options.Noise = false;
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

%Set up a parallel pool (in case it hasn't yet)
nParallelWorkers = feature('numcores');
currPool = gcp('nocreate');
if isempty(currPool)
    try %Handle occasional errors..                                         %"Error using parpool (line 145)
        currPool = parpool('local',nParallelWorkers);                       %Parallel pool failed to start with the following error. For more detailed
    catch %Try again..                                                      %information, validate the profile 'local' in the Cluster Profile Manager.
        pause(5);                                                           %Error using parallel.internal.pool.InteractiveClient>iThrowWithCause (line 670)
        myCluster = parcluster('local');                                    %Failed to initialize the interactive session. Error using  
        delete(myCluster.Jobs);                                             %parallel.internal.pool.InteractiveClient>iThrowIfBadParallelJobStatus (line 810)    
        currPool = parpool('local',nParallelWorkers);                       %The interactive communicating job finished with no message."  
    end                                                                     
end                                                                         %Also: Remove all jobs created with profile local (these are saved after a crash)

%Correct the number of samples such that it is a multiple of nWalkers
samplesPerWalker = ceil(N/K);
mccountRounded = K*samplesPerWalker;  

%Implement parallelization: separate calls to separate cores per round
nCallsPerRound = min(mccountRounded,nParallelWorkers);
samplesPerWalkerPerCall = floor(samplesPerWalker/nCallsPerRound)*ones(1,nCallsPerRound);
samplesPerWalkerPerCall(1) = samplesPerWalkerPerCall(1) + rem(samplesPerWalker,floor(samplesPerWalker/nCallsPerRound));     %The first core/worker may do a few more samples than the others
mccountRoundedPerCall = samplesPerWalkerPerCall*K;

%Covariance Matrix or Variational Posterior available?  --> Careful, we don't check for correctness of the input parameters.. 
covmat = [];
VP = [];
if options.VPbasedStepOut
    VP = widths;
else
    meanX0 = mean(x0,1);                                    %Mean and covariance are used to sample more X0
    if ~isempty(widths) 
        %Obtain covariance matrix from WIDTHS input 
        if numel(widths) == 1; widths = repmat(widths, [1,nDim]); end
        if isvector(widths)
            covmat = diag(widths.^2);                       %Assume no correlations
        else
            covmat = widths;                                %Covariance matrix given by user
        end
    elseif (options.Burnin == 0) && isempty(mcmcStructureLastTime)
        warning('No WIDTHS provided and no BurnIn period requested. But we require a BurnIn to estimate the covmat. So we perform a short N/2 BurnIn anyway');
        options.Burnin = N/2;
    end
end

%Prepare BurnIn round
if options.Burnin > 0
    BurnInsamplesPerWalker = ceil(options.Burnin/K);
    BurninCountRounded = K*BurnInsamplesPerWalker;
    nCallsBurnIn = min(BurninCountRounded,nParallelWorkers);
    BurnInCountPerWalkerPerCall = floor(BurnInsamplesPerWalker/nCallsBurnIn)*ones(1,nCallsBurnIn);
    BurnInCountPerWalkerPerCall(1) = BurnInCountPerWalkerPerCall(1) + rem(BurnInsamplesPerWalker,floor(BurnInsamplesPerWalker/nCallsBurnIn));     %The first core/worker may do a few more samples than the others
    BurninCountRoundedPerCall = BurnInCountPerWalkerPerCall*K;
    BurnInCopy = sum(BurninCountRoundedPerCall);
else
    BurnInCopy = 0;    
end

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
    
    %Overwrite some other settings
    if isempty(VP)
        widths = mcmcStructureLastTime.output.covmat;
        covmat = widths;
    end
    BurnInCopy = mcmcStructureLastTime.output.burn;
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
    funccount = [0 0];
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
        
        %BurnIn round (to estimate the covariance matrix of the posterior)
        if (MCMCroundCounter == 1) && (i_thin == 1) && options.Burnin > 0
            output_BurnIn = cell(1,nCallsBurnIn);                         
            seed_offset = randi(floor(intmax/10));                          %Create a random seed
            parfor i_call=1:nCallsBurnIn
                rng(i_call + seed_offset);                                  %We have to randomize within the parfor loop as the rng seed is otherwise set to "default" on every worker.
                if ~isempty(VP)
                    minit = vbmc_rnd(VP,K,1,1);                             %Sample from variational posterior (VBMC) 
                elseif ~isempty(covmat)
                    minit = mvnrnd(meanX0,covmat,K);                        %Sample from multivariate normal distribution
                else
                    minit = x0;                                             %Don't sample new X0, initialize each BurnIn call with the same samples (not preferable!)
                end
                optionsTmp = options;                                       %Change the options.Burnin setting
                optionsTmp.Burnin = BurninCountRoundedPerCall(i_call);
                %Call EISSAMPLE (by Luigi Acerbi, https://github.com/lacerbi/eissample/). The function was adjusted by David Meijer
                [~,~,~,output_BurnIn{1,i_call}] = eissample_lite_DM2(logPfuns,minit,1,K,widths,LB,UB,optionsTmp,varargin{:});                  %After burn-in period we only collect 1 MCMC sample (which we throw away)!
            end
            %Compute average covariance matrix to set the widths for all following calls of the MCMC algorithm    
            covMats = nan(P.nParams2Fit,P.nParams2Fit,nCallsBurnIn);
            for i_call=1:nCallsBurnIn
                covMats(:,:,i_call) = output_BurnIn{1,i_call}.covmat;       
                funccount = funccount + output_BurnIn{1,i_call}.funccount;                      %Also update the function counts
                nslicecollapsed = nslicecollapsed + output_BurnIn{1,i_call}.nslicecollapsed;
                nTotSamples = nTotSamples + output_BurnIn{1,i_call}.totN;
                CountTemp = sum(~isnan(output_BurnIn{1,i_call}.SigmaFactor));
                SigmaFactor = (SigmaFactorCount*SigmaFactor + CountTemp*nanmean(output_BurnIn{1,i_call}.SigmaFactor)) / (SigmaFactorCount + CountTemp);
                SigmaFactorCount = SigmaFactorCount + CountTemp;
            end
            widths = mean(covMats,3);                                       %The widths are set to the covariance matrix that was estimated as mean covariance matrix across all calls (to parallel workers)       
            covmat = widths;                                                %The covariance matrix is also used for sampling new X0
            options.Burnin = 0;                                             %Reset the burnin options setting (avoid using a Burnin period on consecutive calls of EISSAMPLE)
        end %End of BurnIn round
        
        %Create temporary cell-arrays to store the results
        mcmcParams_Tmp = cell(1,nCallsPerRound);                            %sampled parameters
        mcmcLP_Tmp = cell(1,nCallsPerRound);                                %log-posterior probabilities
        mcmcLPrior_Tmp = cell(1,nCallsPerRound);                            %log-prior probabilities
        mcmcLL_responses_Tmp = cell(1,nCallsPerRound);                      %LL of each sample per set of sampled parameters
        output_Tmp = cell(1,nCallsPerRound);                                %output structure for each call
        
        %Create a new random seed
        seed_offset = randi(floor(intmax/10));
        
        %For each round of sampling we divide the total number of samples into several calls with different initializations  
        parfor i_call=1:nCallsPerRound 
        
            %We have to randomize within the parfor loop as the rng seed is otherwise set to "default" on every worker.
            rng(i_call + seed_offset);
            
            %Sample starting positions for the chains (new on each call)
            if ~isempty(VP)
                minit = vbmc_rnd(VP,K,1,1);                                 %Sample from variational posterior (VBMC) 
            else
                minit = mvnrnd(meanX0,covmat,K);                            %Sample from multivariate normal distribution
            end

            %Call EISSAMPLE (by Luigi Acerbi, https://github.com/lacerbi/eissample/). The function was adjusted by David Meijer
            [mcmcParams_Tmp{1,i_call},mcmcLP_Tmp{1,i_call},~,output_Tmp{1,i_call},mcmcLL_responses_temp] = eissample_lite_DM2(logPfuns,minit,mccountRoundedPerCall(i_call),K,widths,LB,UB,options,varargin{:});
            mcmcLPrior_Tmp{1,i_call} = mcmcLL_responses_temp{1};            
            mcmcLL_responses_Tmp{1,i_call} = mcmcLL_responses_temp{2};             
            
        end %End of parfor loop (CallsPerRound)
        
        %Create an output structure
        if (MCMCroundCounter == 1)
            output = output_Tmp{1,1};
        end
        
        %Update the function counts  
        for i_call=1:nCallsPerRound
            funccount = funccount + output_Tmp{1,i_call}.funccount;                      
            nslicecollapsed = nslicecollapsed + output_Tmp{1,i_call}.nslicecollapsed;
            nTotSamples = nTotSamples + output_Tmp{1,i_call}.totN;
            CountTemp = sum(~isnan(output_Tmp{1,i_call}.SigmaFactor));
            SigmaFactor = (SigmaFactorCount*SigmaFactor + CountTemp*nanmean(output_Tmp{1,i_call}.SigmaFactor)) / (SigmaFactorCount + CountTemp);
            SigmaFactorCount = SigmaFactorCount + CountTemp;
        end
        
        %Concatenate the various calls (also turns cell array into simple matrix)   
        mcmcParams_Tmp = cat(1,mcmcParams_Tmp{:}); 
        mcmcLP_Tmp = cat(1,mcmcLP_Tmp{:}); 
        mcmcLPrior_Tmp = cat(1,mcmcLPrior_Tmp{:}); 
        mcmcLL_responses_Tmp = cat(1,mcmcLL_responses_Tmp{:}); 
        
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
            disp('MCMC did not converge in either of the rounds! Check the data and try again later using different settings (e.g. increase "maxRounds.MCMC").');
        end
        MCMCroundCounter = MCMCroundCounter+1;  %Potential other way to leave the while loop
    end 
end %End of while loop: converged? or maxRounds reached?

%Clear the parallel pool (to get rid of possible lingering memory leaks: https://uk.mathworks.com/matlabcentral/answers/332792-how-to-avoid-memory-leaks-when-function-inside-parfor-generates-warning)
delete(currPool);

%Subtract one counter if we exited the while loop because MCMCroundCounter > P.maxRounds.MCMC (new round was not actually executed)   
if ~MCMCconvergedBool
    MCMCroundCounter = MCMCroundCounter-1;
end

%Check that eissample_lite did run
if ~exist('exitflag','var')
    warning('Eissample did not run! Presumably because "MCMCroundCounter > maxRounds". Increase maxRounds to continue running an old MCMC sequence.');
    [exitflag,R,Neff,tau] = samplediagnose_DM(mcmcParams,K,funccount,nslicecollapsed,trace,options.Noise);  %Ensure that all the necessary output variables exist
end

%Update the output structure
output.nParallelWorkers = nParallelWorkers;         %This is the only additional output parameter in comparison with (non-parallel) eissample_wrapper.m 
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
mcmcStructure.prob.mcmcLogLikelihood = mcmcLP-mcmcLPrior;                    %Should be equal to "sum(mcmcLL_responses,2)"
mcmcStructure.prob.mcmcLogPosterior = mcmcLP;
mcmcStructure.prob.mcmcLL_responses = mcmcLL_responses;

mcmcStructure.converge.exitflag = exitflag;
mcmcStructure.converge.maxR = max(R);
mcmcStructure.converge.minNeff = min(Neff);

mcmcStructure.output = output;

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


