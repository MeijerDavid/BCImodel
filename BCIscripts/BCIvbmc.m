function dataFit_VBMC = BCIvbmc(BCIfitResults)
%Perform Variational Bayesian Monte Carlo to obatin a rough estimate of the
%posterior probability distribution over parameters.

%Start a timer
disp('Starting VBMC fit ...');
cStart = clock;  

%Add relevant 'VBMC' folders to the path
me = mfilename;                                             %what is my filename
pathstr = fileparts(fileparts(which(me)));                  %get location of 'BCImodel' folder
addpath(genpath([pathstr filesep 'vbmc']));  

%Get the data and some settings
[P,trueLocsAV,responsesAVC,i_Conditions] = unpackBCIfitResults(BCIfitResults);
LB = P.Bounds.conv.LB;
UB = P.Bounds.conv.UB;
PLB = P.Bounds.conv.PLB;
PUB = P.Bounds.conv.PUB;

%Get the optimized parameters that were fitted by BADS
params0 = BCIconvertParams2Fit(BCIfitResults.BADS.fittedParams,P,'real2conv');

%Create anonymous fitting function for input to VBMC
LLfun = @(params) BCIcompLL(params,P,trueLocsAV,responsesAVC,i_Conditions);               %Returns a log-likelihood vector (one LL value for each response)
Probfun = @(params) sum(LLfun(params))+BCIcomputeLogPrior(params,P);                      %Returns the total posterior probability (LL + logPrior)

%Options for VBMC
VBMCoptions = vbmc('defaults');
VBMCoptions.Display = 'off';                                        %No progress in Command Window please
VBMCoptions.MaxIter = 100*(2+P.nParams2Fit);                        %Max number of iterations (default is 50*(2+nParams));
VBMCoptions.MaxFunEvals = 100*(2+P.nParams2Fit);                    %Max number of target fcn evaluations (default is 50*(2+nParams))';
VBMCoptions.RetryMaxFunEvals = 100*(2+P.nParams2Fit);               %If not converged, VBMC will internally retry with a new initialization based on the non-converged solution

%Call VBMC wrapper function
if P.parallel.VBMC
    dataFit_VBMC = vbmc_wrapper_par(Probfun,params0,LB,UB,PLB,PUB,VBMCoptions,P.nAttempts.VBMC);        %Make use of Parallel Computing Toolbox
else
    dataFit_VBMC = vbmc_wrapper(Probfun,params0,LB,UB,PLB,PUB,VBMCoptions,P.nAttempts.VBMC);            %Default
end

%Report computation time in command window
fprintf('Finished VBMC fit, elapsed time (days hours:minutes:seconds) %s \n',datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

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
