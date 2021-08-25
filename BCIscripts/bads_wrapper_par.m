function badsStructure = bads_wrapper_par(fun,x0,LB,UB,PLB,PUB,nonbcon,options,nAttempts,nGrid,varargin)
%Small wrapper around BADS in which we call BADS multiple times and attempt
%to select the overall minimum rather than a local minimum. We make use of
%parallel computing. A gridsearch is performed before BADS is called.

%Initialize output
badsStructure = [];

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

%Find number of parameters to fit
equalNparams = (length(LB) == length(UB)) && (length(LB) == length(PLB)) && (length(LB) == length(PUB));
assert(equalNparams,'all bounds need to have equal length'); 
nParams = length(LB);
if isempty(x0)
    x0 = nan(1,nParams);
else
    assert(length(x0) == nParams,'x0 needs to have same length as bounds, or set as empty');
end

%Set defaults
if nargin < 9
    nAttempts = [1 4];
end
if nargin < 10
    nGrid = nParams*1000;
end

%Assure that there are enough gridpoints to start the required number of tasks   
nTasks2Start = nAttempts(2)+nParallelWorkers;
nGrid = max(nGrid,nTasks2Start);

%Pass varargin into objective function
if isempty(varargin)
    funwrapper = fun;   % No additional function arguments passed
else
    funwrapper = @(x) fun(x,varargin{:});
end

%Randomly draw a 'grid' for the gridsearch
params_grid = nan(nGrid,nParams);
for i=1:nParams
    if isnan(x0(i))
        params_grid(:,i) = (PUB(i)-PLB(i))*rand([nGrid,1])+PLB(i);                  %Randomly draw from a uniform distribution between the plausible bounds
    else
        params_grid(:,i) = sampleXnormal(nGrid,x0(i),PLB(i),PUB(i),LB(i),UB(i));    %Randomly draw from a normal distribution centred on x0 and with SD relative to the plausible bounds
    end
end

%Perform grid search on the specified parameter space
tic;
LP_gridSearch = nan(nGrid,1);
parfor k = 1:nGrid
    LP_gridSearch(k) = funwrapper(params_grid(k,:));
end
disp(['Gridsearch (nGrid = ' num2str(nGrid) ') finished in: ' datestr(toc/(24*60*60),'HH:MM:SS') ' (HH:MM:SS).']);

%Choose starting positions from the grid search
[~,idx_sorted] = sort(LP_gridSearch);                                                               %Ascending order is default (i.e. lowest first. Function values will be minimized by BADS, so smallest values are "best")
params0_tmp = params_grid(idx_sorted(1:nTasks2Start),:);

%Call BADS as many times as maximally needed to keep all workers busy all the time 
FuturesArray(nTasks2Start,1) = parallel.FevalFuture();
for k = 1:nTasks2Start
    FuturesArray(k) = parfeval(@bads,4,funwrapper,params0_tmp(k,:),LB,UB,PLB,PUB,nonbcon,options);  %Expect 4 output arguments, and use 8 input arguments
end
idx_running = 1:nTasks2Start;

%Initialize BADS outputs
params0 = nan(nAttempts(2),nParams);
fittedParams = nan(nAttempts(2),nParams);
fittedValue = nan(nAttempts(2),1);
exitflags = nan(nAttempts(2),1);
outputstructs = cell(nAttempts(2),1);

%Attempt to optizime the parameters, repeat minimally nAttempts(1) times and, if necessary, maximally nAttempts(2) times   
BADSconvergedBool = false;
nFinished = 0;
while ~BADSconvergedBool && (nFinished < nAttempts(2))
    
    %Update the counter (prematurely)
    nFinished = nFinished+1;
    
    %Retrieve the BADS results from the parallel workers as they become available    
    [completedIdx,fittedParams(nFinished,:),fittedValue(nFinished),exitflags(nFinished),outputstructs{nFinished,1}] = fetchNext(FuturesArray); 
    params0(nFinished,:) = params0_tmp(idx_running(completedIdx),:);
    FuturesArray(completedIdx) = [];    %Release memory (https://uk.mathworks.com/matlabcentral/answers/424473-parfeval-memory-consumption-piling-up-clear-output-data)
    idx_running(completedIdx) = [];
    
    %After enough BADS searches have been performed...
    if any(nFinished >= nAttempts)
        
        %Choose final best fit (highest LL = lowest negLL = lowest fittedValue)
        [~,idx_sorted] = sort(fittedValue(1:nFinished,1)); %Ascending order..
        BestIdxCounter = 1;
        while ~BADSconvergedBool && (BestIdxCounter <= nFinished)
            idxTemp = idx_sorted(BestIdxCounter);
            if exitflags(idxTemp,1) > 0
                bestIdx = idxTemp;
                BADSconvergedBool = true;
            end
            BestIdxCounter = BestIdxCounter+1;
        end

        %If maximum attempts reached and not converged, then best result is the non-converged lowest negLL  
        if ~BADSconvergedBool && (nFinished >= nAttempts(2))
            disp(['BADS fit did not converge in either of the ' num2str(nFinished) ' runs! Check the data and try again later using different settings...']);
            bestIdx = idx_sorted(1);
            break; %Break from while loop
        end
    end
    
end %End of while loop (BADS converged?)

%Cancel any remaining jobs on the cores (running or queued)
cancel(FuturesArray);
delete(FuturesArray);       %Delete the objects
clear FuturesArray          %Clear the references to those objects (https://uk.mathworks.com/help/parallel-computing/parallel.job.delete.html)    

%Clear the parallel pool (to get rid of possible lingering memory leaks: https://uk.mathworks.com/matlabcentral/answers/332792-how-to-avoid-memory-leaks-when-function-inside-parfor-generates-warning)
delete(currPool);

%Collect best results to return in struct
badsStructure.wrapper.nParallelWorkers = nParallelWorkers;
badsStructure.wrapper.nFinished = nFinished;
badsStructure.wrapper.bestIdx = bestIdx;
badsStructure.wrapper.params0 = params0(1:nFinished,:);
badsStructure.wrapper.fittedParams = fittedParams(1:nFinished,:);
badsStructure.wrapper.fittedValue = fittedValue(1:nFinished,1);
badsStructure.wrapper.exitflags = exitflags(1:nFinished,1);

badsStructure.fittedParams = fittedParams(bestIdx,:);
badsStructure.fittedValue = fittedValue(bestIdx);
badsStructure.exitflag = exitflags(bestIdx);
badsStructure.output = outputstructs{bestIdx,1};

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function X = sampleXnormal(N,meanX,PLB,PUB,LB,UB)
%Sample from univariate normal distribution with SD set as: [PUB-PLB]/6. 
%Ensure that all samples fall within LB and UB

%Check that meanX is within the bounds
if (meanX < LB) || (meanX > UB)
    error('Error when creating intializations for BADS: meanX falls outside of bounds [LB UB]');
end

%Set Standard Deviation
SD = (PUB-PLB)/6;

%Sample until all X fall within bounds
X = nan(N,1);
invalid = true(N,1);
maxAttempts = 1e5;
attemptCounter = 0;
while any(invalid) && (attemptCounter < maxAttempts)
    attemptCounter = attemptCounter+1;
    idx_invalid = find(invalid);
    X(idx_invalid,1) = normrnd(meanX,SD,[numel(idx_invalid) 1]);
    invalid = (X <= LB) | (X >= UB);
end

%Throw error if not all valid
if any(invalid)
    error('Error when creating intializations for BADS: unable to find X within bounds [LB UB]');
end

end %[EoF] 
