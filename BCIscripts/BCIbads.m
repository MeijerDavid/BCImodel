function dataFit_BADS = BCIbads(BCIfitResults)
%Perform BADS parameter optimization

%Start a timer
disp('Starting BADS fit ...');
cStart = clock;  

%Add relevant 'BADS' folders to the path
me = mfilename;                                      %what is my filename
pathstr = fileparts(fileparts(which(me)));           %get location of 'BCImodel' folder
addpath(genpath([pathstr filesep 'bads']));          %add 'BADS' folder to the path
addpath(genpath([pathstr filesep 'gofit']));         %add 'Goodness-of-Fit' folder to the path

%Get the data and some settings
[P,trueLocsAV,responsesAVC,i_Conditions] = unpackBCIfitResults(BCIfitResults);
LB = P.Bounds.conv.LB;
UB = P.Bounds.conv.UB;
PLB = P.Bounds.conv.PLB;
PUB = P.Bounds.conv.PUB;

%Create anonymous function for input to BADS
LLfun = @(params) BCIcompLL(params,P,trueLocsAV,responsesAVC,i_Conditions);                         %Returns a log-likelihood vector (one LL value for each response)
if P.ForceLogLikelihood
    Probfun = @(params) sum(LLfun(params));                                                         %Returns the total log-likelihood 
else
    Probfun = @(params) sum(LLfun(params))+BCIcomputeLogPrior(params,P);                            %Returns the total posterior probability (LL + logPrior)
end
fitfun = @(params) -Probfun(params);                                                                %Returns the negative log probability

%Set initialization parameters for grid search
params0 = setInitialization(P);

%Set some options for BADS
BADSoptions = bads('defaults');
BADSoptions.Display = 'off';

%Call BADS wrapper function
if P.parallel.BADS
    dataFit_BADS = bads_wrapper_par(fitfun,params0,LB,UB,PLB,PUB,[],BADSoptions,P.nAttempts.BADS,P.nGridSearch);        %Make use of Parallel Computing Toolbox
else
    dataFit_BADS = bads_wrapper(fitfun,params0,LB,UB,PLB,PUB,[],BADSoptions,P.nAttempts.BADS,P.nGridSearch);            %Default
end 

%Collect the probabilities (prior, likelihood, posterior)
fittedValue = -dataFit_BADS.fittedValue;                                    %Note change from negative log(prob) to positive log(prob)
dataFit_BADS = rmfield(dataFit_BADS,'fittedValue');
prob.logPrior = BCIcomputeLogPrior(dataFit_BADS.fittedParams,P);
if P.ForceLogLikelihood
    prob.logLikelihood = fittedValue;
    prob.logPosterior = fittedValue + prob.logPrior;
else
    prob.logLikelihood = fittedValue - prob.logPrior;
    prob.logPosterior = fittedValue;
end
dataFit_BADS.prob = prob;

%Compute Goodness-of-Fit (use original settings and conditions structure - not the collapsed one, see "unpackBCIfitResults")  
dataFit_BADS.goodnessOfFit = BCIcomputeGoF(BCIfitResults.settings,trueLocsAV,responsesAVC,BCIfitResults.data.i_Conditions,'BADS',dataFit_BADS.fittedParams);

%Convert the fitted parameters to normal space
dataFit_BADS.fittedParams = BCIconvertParams2Fit(dataFit_BADS.fittedParams,P,'conv2real');
    
%Report computation time in command window
fprintf('Finished BADS fit, elapsed time (days hours:minutes:seconds) %s \n',datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end %[EoF]

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
    binIdxResponsesAVC((responsesAVC(:,3) == 1),3) = 1;                                 %A "1" means common-source
    binIdxResponsesAVC((responsesAVC(:,3) ~= 1) & ~isnan(responsesAVC(:,3)),3) = 2;     %Any other (e.g. "0" or "2") but not NaN means independent sources!
    
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

function params0 = setInitialization(P)

    %Initialize the parameters with the user defined values (default values will have been overwritten with NaNs for parameters that are fitted)    
    params0 = nan(1,P.nParams2Fit);
    for i=1:P.nParams2Fit
        params0(i) = P.(P.ParamNames2Fit{1,i});
    end

    %Convert the parameters (to log or logit space, in which we fit the parameters)
    params0 = BCIconvertParams2Fit(params0,P,'real2conv');

    %Ensure that the requested parameter initializations are within the bounds!
    for i=1:P.nParams2Fit
        if params0(i) <= P.Bounds.conv.LB(i)
            error(['Requested parameter initialization ' num2str(i) ' is smaller/equal to hard lower bound']);
        elseif params0(i) >= P.Bounds.conv.UB(i)
            error(['Requested parameter initialization ' num2str(i) ' is larger/equal to hard upper bound']);
        elseif params0(i) < P.Bounds.conv.PLB(i)
            warning(['Requested parameter initialization ' num2str(i) ' is smaller than plausible lower bound. It will be moved to the plausible lower bound...']);
            params0(i) = P.Bounds.conv.PLB(i);
        elseif params0(i) > P.Bounds.conv.PUB(i)
            warning(['Requested parameter initialization ' num2str(i) ' is larger than plausible upper bound. It will be moved to the plausible upper bound...']);
            params0(i) = P.Bounds.conv.PUB(i);
        end
    end
    
end %[EoF]
