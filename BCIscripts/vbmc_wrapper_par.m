function vbmcStructure = vbmc_wrapper_par(logPfun,x0,LB,UB,PLB,PUB,options,nAttempts,varargin)
%Small wrapper around VBMC in which we call VBMC multiple times and attempt
%to reach convergence. We will also validate the result across the various
%VPs that are obtained. We make use of parallel computing. 

%Initialize output
vbmcStructure = [];

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

%Set default number of rounds
if nargin < 8
    nAttempts = [1 4];
end

%Pass varargin into objective function
if isempty(varargin)
    logPfunwrapper = logPfun;   % No additional function arguments passed
else
    logPfunwrapper = @(x) logPfun(x,varargin{:});
end

%Prepare initializations for the VBMC procedure  
nTasks2Start = nAttempts(2)+nParallelWorkers;
jitterFactorSmaller = 100;
[meanX,jitterPLB,jitterPUB] = shrinkPlausibleBounds(x0,PLB,PUB,LB,UB,jitterFactorSmaller);
init_tmp = sampleXnormal(nTasks2Start,meanX,jitterPLB,jitterPUB);                   %Sample initializations from multivariate normal around x0
vbmcFactorSmaller = 10;
[~,vbmcPLB,vbmcPUB] = shrinkPlausibleBounds(x0,PLB,PUB,LB,UB,vbmcFactorSmaller);    %The smaller interval for the plausible bounds helps the VBMC algorithm and avoids a warning when init is outside of the plausible bounds
                                                    
%Call VBMC as many times as maximally needed to keep all workers busy all the time 
FuturesArray(nTasks2Start,1) = parallel.FevalFuture();
for k = 1:nTasks2Start
    FuturesArray(k) = parfeval(@vbmc,5,logPfunwrapper,init_tmp(k,:),LB,UB,vbmcPLB,vbmcPUB,options);  %Expect 5 output arguments, and use 8 input arguments
end
idx_running = 1:nTasks2Start;

%Initialize arrays for the finished results
init = nan(size(init_tmp));
vp_all = cell(nAttempts(2),1);
elbo_all = nan(nAttempts(2),1);
elbo_sd_all = nan(nAttempts(2),1);
exitflag_all = nan(nAttempts(2),1);
output_all = cell(nAttempts(2),1);

%Initialize counters
i_error = false(nTasks2Start,1);
nErrors = 0;

nFinished = 0;
VBMCconvergedBool = 0;

%Try to find a coverged solution; repeat minimally nAttempts(1) times and, if necessary, maximally nAttempts(2) times   
while ~VBMCconvergedBool && (nFinished < nAttempts(2))
    
    %Update the counter (prematurely)
    nFinished = nFinished+1;
            
    %Retrieve the VBMC results from the parallel workers as they become available    
    try 
        [completedIdx,vp_all{nFinished,1},elbo_all(nFinished,1),elbo_sd_all(nFinished,1),exitflag_all(nFinished,1),output_all{nFinished,1}] = fetchNext(FuturesArray); 
        init(nFinished,:) = init_tmp(idx_running(completedIdx),:);
    catch ME                                                                %We run fetchNext in a try-catch because VBMC v1.0 occasionally resulted in an error (probably an out-of-bounds bug somewhere..)
        completedIdx = find([FuturesArray.Read],1);
        i_error(nFinished) = true;
        nErrors = nErrors+1;
        initErrors(nErrors,:) = init_tmp(idx_running(completedIdx),:);
        msgStructErrors{nErrors} = ME;
        warning(['One of the VBMC attempts returned with an error. The error message will be saved. Init: ' num2str(initErrors(nErrors,:))]);
    end
    FuturesArray(completedIdx) = [];    %Release memory (https://uk.mathworks.com/matlabcentral/answers/424473-parfeval-memory-consumption-piling-up-clear-output-data)
    idx_running(completedIdx) = [];
    
    %After enough VBMC solutions have been obtained...
    if any(nFinished >= nAttempts)
        
        %Check whether any one of the VBMC fits converged
        VBMCconvergedBool = any(exitflag_all);
        
        %If maximum attempts reached and not converged, then break from the while loop  
        if ~VBMCconvergedBool && (nFinished >= nAttempts(2))
            disp(['VBMC fit did not converge in either of the ' num2str(nFinished) ' runs! Check the data and try again later using different settings...']);
            break; %break from while loop
        end
    end
 
end %End of while loop (VBMC converged?)

%Cancel any remaining jobs on the cores (running or queued)
cancel(FuturesArray);
delete(FuturesArray);       %Delete the objects
clear FuturesArray          %Clear the references to those objects (https://uk.mathworks.com/help/parallel-computing/parallel.job.delete.html)    

%Clear the parallel pool (to get rid of possible lingering memory leaks: https://uk.mathworks.com/matlabcentral/answers/332792-how-to-avoid-memory-leaks-when-function-inside-parfor-generates-warning)
delete(currPool);

%Deal with errors
if nErrors > 0
    init(i_error,:) = [];
    vp_all(i_error,:) = [];             %Note that i_error may be longer than vp_all, but this is not a problem: length(i_error) = nTasks2Start, length(vp_all) = nAttempts(2)
    elbo_all(i_error,:) = [];
    elbo_sd_all(i_error,:) = [];
    exitflag_all(i_error,:) = [];
    output_all(i_error,:) = [];
    nFinished = nFinished - nErrors;
end

if nFinished == 0 
    %Create an output struct with information about the errors
    warning('All VBMC attempts resulted in an error. Check the data and try again.');
    vbmcStructure.nParallelWorkers = nParallelWorkers;
    vbmcStructure.nFinished = nFinished;
    vbmcStructure.nErrors = nErrors;
    vbmcStructure.params0 = initErrors;
    vbmcStructure.msgStructErrors = msgStructErrors;

else
    %Validate the results 
    [~,~,idx_best,stats_Validation] = vbmc_diagnostics_DM(vp_all(1:nFinished,1));

    %Find the mean and SDs for each parameter estimate of all VPs
    [means,SDs] = VPparamSummary(vp_all(1:nFinished,1));

    %Save in FitResults structure for output
    vbmcStructure.wrapper.nParallelWorkers = nParallelWorkers;
    vbmcStructure.wrapper.nFinished = nFinished;
    vbmcStructure.wrapper.idx_best = idx_best;
    vbmcStructure.wrapper.params0 = init(1:nFinished,:);
    vbmcStructure.wrapper.paramMeans = means;
    vbmcStructure.wrapper.paramSDs = SDs;
    vbmcStructure.wrapper.vp = vp_all(1:nFinished,1);
    vbmcStructure.wrapper.elbo = elbo_all(1:nFinished,1);
    vbmcStructure.wrapper.elbo_sd = elbo_sd_all(1:nFinished,1);
    vbmcStructure.wrapper.exitflag = exitflag_all(1:nFinished,1);

    if nErrors > 0
        vbmcStructure.errors.nErrors = nErrors;
        vbmcStructure.errors.params0 = initErrors;
        vbmcStructure.errors.msgStructErrors = msgStructErrors;
    end
    
    if ~isempty(idx_best)                                                           %idx_best will be empty if none of the VBMC attempts converged
        vbmcStructure.vp = vp_all{idx_best,1};
        vbmcStructure.elbo = elbo_all(idx_best,1);
        vbmcStructure.elbo_sd = elbo_sd_all(idx_best,1);
        vbmcStructure.exitflag = exitflag_all(idx_best,1);
        vbmcStructure.output = output_all{idx_best,1};

        vbmcStructure.validation = stats_Validation;

        %Finally, compute some stats for the best VP
        Xs = vbmc_rnd(vbmcStructure.vp,1e4,1,1);                                    %Generate ten thousand samples from the variational posterior (conservative with memory, you could choose more samples..)    
        vbmcStructure.stats.quantiles = [0.025,0.05,0.25,0.5,0.75,0.95,0.975];      %We'll ask for 7 quantiles to describe the distribution of each parameter
        vbmcStructure.stats.post_iqr = quantile(Xs,vbmcStructure.stats.quantiles,1);%posterior quantiles (0.5 = median). This returns a matrix with size [7 x nParams]
        vbmcStructure.stats.post_mean = mean(Xs,1);                                 %posterior mean for each parameter
        vbmcStructure.stats.post_std = std(Xs,[],1);                                %posterior standard deviation for each parameter
        vbmcStructure.stats.post_cov = cov(Xs);                                     %posterior covariance matrix in real space

        % %You can visualize the posterior marginals using the CORNERPLOT function
        % cornerplot(Xs);
    end
end

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function [meanX,newPLB,newPUB] = shrinkPlausibleBounds(meanX,PLB,PUB,LB,UB,factorSmaller)
%Decrease the space within the plausible bounds by some factorSmaller 
%Use meanX as centre point for the new plausible bounds. 
%Also ensure that meanX is within the hard bounds.

%Size of current/old plausible bounds
rangeOld = PUB-PLB;
rangeNew = rangeOld/factorSmaller;
rangeEitherSide = rangeNew/2;

%Check that meanX is within the hard bounds (i.e. not ON the hard bounds),
%if on/outside, move meanX to halfway the hard and current/old plausible bounds    
if any(meanX <= LB) || any(meanX >= UB)
    warning('x0 is outside or on bounds LB/UB. x0 will be moved to halfway hard and plausible bounds.');
end
meanX(meanX <= LB) = mean([LB(meanX <= LB); PLB(meanX <= LB)],1);           
meanX(meanX >= UB) = mean([LB(meanX >= UB); PLB(meanX >= UB)],1);           %Note that we assumed that [meanX,PLB,PUB,LB,UB] are all row vectors     

%Check that meanX is within the current/old plausible bounds, 
%if outside, move the current/old plausible bounds to halfway meanX and the hard bounds    
%Note that these new plausible bounds are only used to sample initializations for the VBMC algorithm. They do not change the trapezoid prior etc..    
PLB(meanX < PLB) = mean([meanX(meanX < PLB); LB(meanX < PLB)],1);           
PUB(meanX > PUB) = mean([meanX(meanX > PUB); UB(meanX > PUB)],1);           %Note that we assumed that [meanX,PLB,PUB,LB,UB] are all row vectors     

%Suggest new plausible bounds (we can only make the plausible space smaller, not larger)
newPLB = max([PLB; meanX-rangeEitherSide],[],1);
newPUB = min([PUB; meanX+rangeEitherSide],[],1);                            %Note that we assumed that [meanX,PLB,PUB] are all row vectors 

end %[EOF]

%-------------------------------------------------------------------------

function X = sampleXnormal(N,meanX,PLB,PUB)
%Sample from multivariate normal distribution with SDs set as: [PUB-PLB]/6. 
%Ensure that all samples fall within PUB and PLB

%Check that meanX is within the bounds
if any(meanX < PLB) || any(meanX > PUB)
    error('Error when creating intializations for VBMC: meanX falls outside of bounds [PLB PUB]');
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
    error('Error when creating intializations for VBMC: unable to find X within bounds [PLB PUB]');
end

end %[EoF] 


%-------------------------------------------------------------------------

%Based on "vbmc_diagnostics" but simplified and without command window output
function [exitflag,best,idx_best,stats] = vbmc_diagnostics_DM(vp_array,beta_lcb,elbo_thresh,sKL_thresh,maxmtv_thresh)
%VBMC_DIAGNOSTICS Convergence diagnostics between multiple VBMC runs.
%   EXITFLAG = VBMC_DIAGNOSTICS(VP_ARRAY) runs a series of diagnostic tests 
%   on an array of variational posteriors. VP_ARRAY is a cell array or 
%   struct array of variational posteriors obtained by separate runs of 
%   VBMC on the same problem. EXITFLAG describes the result of the analysis.
%   Possible values of EXITFLAG and the corresponding test results are
%
%    1  PASSED: All diagnostics tests passed.
%    0  FAILED: Only one solution converged, cannot perform useful diagnostics.
%   -1  FAILED: Not enough solutions agree with the best posterior (in terms 
%       of symmetrized KL-divergence or maximum marginal total variation 
%       distance).
%   -2  FAILED: Not enough solutions agree with the best ELBO.
%   -3  FAILED: No solution has converged. No "best" solution.
%
%   A minimum of 2 separate runs of VBMC are required to perform diagnostic
%   checks, and it is recommended to perform at least 3 or 4 runs.
%
%   [EXITFLAG,BEST] = VBMC_DIAGNOSTICS(...) returns a struct BEST that 
%   contains the "best" solution, that is the solution with highest ELCBO 
%   (lower confidence bound on the ELBO) among the solutions that have 
%   converged. The fields of BEST are 'vp', which contains the variational 
%   posterior, its associated ELBO ('elbo') and estimated error ('elbo_sd'). 
%   You should be wary of using a solution which has not fully passed the 
%   test diagnostics. BEST is returned empty if no solution has converged.
%
%   [EXITFLAG,BEST,IDX_BEST] = VBMC_DIAGNOSTICS(...) also returns the index
%   within the array VP_ARRAY of the returned "best" solution. IDX_BEST
%   is returned empty if no solution has converged.
%
%   [EXITFLAG,BEST,IDX_BEST,STATS] = VBMC_DIAGNOSTICS(...) returns a struct
%   STATS with summary statistics of the diagnostic tests.
%
%   [...] = VBMC_DIAGNOSTICS(VP_ARRAY,BETA_LCB) uses lower confidence bound
%   factor BETA_LCB to judge the best solution in terms of ELCBO (default
%   BETA_LCB = 3).
%
%   [...] = VBMC_DIAGNOSTICS(VP_ARRAY,BETA_LCB,ELBO_THRESH) specifies the
%   threshold on the ELBO difference to judge two variational solutions as 
%   "close" (default ELBO_THRESH = 1).
%
%   [...] = VBMC_DIAGNOSTICS(VP_ARRAY,BETA_LCB,ELBO_THRESH,SKL_THRESH) 
%   specifies the threshold on the symmetrized Kullback-Leibler divergence
%   to judge two variational posteriors as "close" (default SKL_THRESH = 1).
%
%   [...] = VBMC_DIAGNOSTICS(VP_ARRAY,BETA_LCB,ELBO_THRESH,SKL_THRESH,MAXMTV_THRESH) 
%   specifies the threshold on the maximum marginal total variation distance
%   to judge two variational posteriors as "close" (default MAXMTV_THRESH = 0.2).
%
%   See also VBMC, VBMC_EXAMPLES, VBMC_KLDIV, VBMC_MTV.

if nargin < 3 || isempty(beta_lcb); beta_lcb = 3; end
if nargin < 4 || isempty(elbo_thresh); elbo_thresh = 1; end                   
if nargin < 5 || isempty(sKL_thresh); sKL_thresh = 1; end                   
if nargin < 6 || isempty(maxmtv_thresh); maxmtv_thresh = 0.2; end

Nruns = numel(vp_array);
exitflag = Inf;

% At least one third of solutions need to be close to the best
TolClose = 1/3;

% Get stats for each run
elbo = NaN(1,Nruns);    elbo_sd = NaN(1,Nruns);     stable_flag = false(1,Nruns);
for iFit = 1:Nruns
    elbo(iFit) = vp_array{iFit}.stats.elbo;
    elbo_sd(iFit) = vp_array{iFit}.stats.elbo_sd;
    stable_flag(iFit) = vp_array{iFit}.stats.stable;
end

% Check which runs have converged
idx_ok = stable_flag;
idx_active = idx_ok;

% None has converged
if sum(idx_ok) == 0
    idx_active = true(size(idx_ok));       %Use the non-converged solutions for defining which is the best (based on elcbo)
    exitflag = -3;  
% Only one has converged   
elseif sum(idx_ok) == 1
    exitflag = 0;
end

% Compute ELCBO, that is lower confidence bound on ELBO
elcbo = elbo - beta_lcb*elbo_sd;

% Pick best variational solution based on ELCBO
elcbo_eff = elcbo;
elcbo_eff(~idx_active) = -Inf;
[~,idx_best] = max(elcbo_eff);

% Compute distances (KL-divergence and MaxMTV) across all pairs of solutions
kl_mat = zeros(Nruns,Nruns);
maxmtv_mat = zeros(Nruns,Nruns);
for iRun = 1:Nruns
    for jRun = iRun+1:Nruns        
        [kl,xx1,xx2] = vbmc_kldiv(vp_array{iRun},vp_array{jRun});
        kl_mat(iRun,jRun) = kl(1);
        kl_mat(jRun,iRun) = kl(2);
        maxmtv_mat(iRun,jRun) = max(vbmc_mtv(xx1,xx2));         %maximum across dimensions..
        maxmtv_mat(jRun,iRun) = maxmtv_mat(iRun,jRun);
    end
end

% Compute symmetrized KL-divergence between best solution and the others
sKL_best = NaN(1,Nruns);
for iRun = 1:Nruns
    sKL_best(iRun) = 0.5*(kl_mat(iRun,idx_best)+kl_mat(idx_best,iRun));
end

% Max marginal total variation between best solution and the others
maxmtv_best = maxmtv_mat(idx_best,:);

if Nruns > 1
    % Check closeness of solutions in terms of ELBO
    elbo_ok = abs(elbo(idx_best) - elbo) < elbo_thresh;
    if sum(elbo_ok) < max(Nruns*TolClose,2)
       exitflag = min(exitflag,-2);     %Not enough solutions agree with the best ELBO
    end

    % Check closeness of solutions in terms of symmetrized KL-divergence
    sKL_ok = sKL_best < sKL_thresh;
    if sum(sKL_ok) < max(Nruns*TolClose,2)
       exitflag = min(exitflag,-1);     %Not enough solutions agree with the best posterior in terms of symmetrized KL-divergence
    end
    
    % Check closeness of solutions in terms of max MTV
    maxmtv_ok = maxmtv_best < maxmtv_thresh;
    if sum(maxmtv_ok) < max(Nruns*TolClose,2)
       exitflag = min(exitflag,-1);     %Not enough solutions agree with the best posterior in terms of max marginal total variation distance
    end
    
end

% Nothing bad found, diagnostic test passed
if isinf(exitflag); exitflag = 1; end

switch exitflag
    case 1; msg = 'Diagnostic test PASSED.';
    case 0; msg = 'Diagnostic test FAILED. Only one solution converged; cannot perform useful diagnostics.';
    case -1; msg = 'Diagnostic test FAILED. Not enough solutions agree with the best posterior (in terms of symmetrized KL-divergence or maximum marginal total variation distance).';
    case -2; msg = 'Diagnostic test FAILED. Not enough solutions agree with the best ELBO.';
    case -3; msg = 'Diagnostic test FAILED. No solution has converged.';    
end


% Return best solution, only if it has converged
if nargout > 1
    if stable_flag(idx_best)
        best.vp = vp_array{idx_best};    
        best.elbo = elbo(idx_best);
        best.elbo_sd = elbo_sd(idx_best);
    else
        best = [];
        idx_best = [];
    end
end

% Create diagnostics STATS struct
if nargout > 3
    stats.beta_lcb = beta_lcb;
    stats.elbo_thresh = elbo_thresh;
    stats.sKL_thresh = sKL_thresh;
    stats.maxmtv_thresh = maxmtv_thresh;
    stats.elbo = elbo;
    stats.elbo_sd = elbo_sd;
    stats.elcbo = elcbo;
    stats.idx_best = idx_best;
    stats.sKL_best = sKL_best;
    stats.maxmtv_best = maxmtv_best;
    stats.kl_mat = kl_mat;
    stats.maxmtv_mat = maxmtv_mat;
    stats.exitflag = exitflag;
    stats.msg = msg;
end

end %[EoF]

%-------------------------------------------------------------------------

function [means,SDs] = VPparamSummary(vp_all)

nVPs = length(vp_all);
nParams = vp_all{1}.D;

means = nan(nVPs,nParams);
SDs = nan(nVPs,nParams);
for i=1:nVPs
    Xs = vbmc_rnd(vp_all{i},1e4,1,1);                                       %Generate ten thousand samples from the variational posterior    
    means(i,:) = mean(Xs,1);                                                %posterior mean for each parameter
    SDs(i,:) = std(Xs,[],1);                                                %posterior standard deviation for each parameter
end

end %[EoF]

