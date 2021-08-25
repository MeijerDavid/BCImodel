function badsStructure = bads_wrapper(fun,x0,LB,UB,PLB,PUB,nonbcon,options,nAttempts,nGrid,varargin)
%Small wrapper around BADS in which we call BADS multiple times and attempt
%to select the overall minimum rather than a local minimum. A gridsearch is
%performed before BADS is called.

%Initialize output
badsStructure = [];         

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
nGrid = max(nGrid,nAttempts(2));

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
%profile on
tic;
LP_gridSearch = nan(nGrid,1);
%T_gridSearch = nan(nGrid,1);
for k = 1:nGrid
    %tic;
    %profile on
    %for i=1:100
        LP_gridSearch(k) = funwrapper(params_grid(k,:));
    %end
    %profile viewer
    %T_gridSearch(k) = toc;
end
disp(['Gridsearch (nGrid = ' num2str(nGrid) ') finished in: ' datestr(toc/(24*60*60),'HH:MM:SS') ' (HH:MM:SS).']);
%profile viewer

%Choose starting positions from the grid search
[~,idx_sorted] = sort(LP_gridSearch);                                               %Ascending order is default (i.e. lowest first. Function values will be minimized by BADS, so smallest values are "best")
params0 = params_grid(idx_sorted(1:nAttempts(2)),:);

%Initialize BADS outputs
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
    
    %Call BADS minimization algorithm
    [fittedParams(nFinished,:),fittedValue(nFinished),exitflags(nFinished),outputstructs{nFinished,1}] = bads(funwrapper,params0(nFinished,:),LB,UB,PLB,PUB,nonbcon,options);
    
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

%Collect best results to return in struct
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
