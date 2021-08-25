function P = BCIsetOptions(OptionsStruct,responsesAVC)
% Set defaults and fixed values for parameter settings and parameter 
% bounds for the fitting algorithm. Also see "BCIsetDefaults.m"

P = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set the parameter values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the names of the parameters to fit
if isfield(OptionsStruct,'ParamNames2Fit')
    P.ParamNames2Fit = OptionsStruct.ParamNames2Fit;
    OptionsStruct = rmfield(OptionsStruct,'ParamNames2Fit');
else
    P.ParamNames2Fit = {'Pcommon','sigmaP','sigmaA','sigmaV'};              %default: 4 standard params of BCI model
end
P.nParams2Fit = numel(P.ParamNames2Fit);

%Get the condition numbers of the parameters to fit
if isfield(OptionsStruct,'ParamsPerCond')
    P.ParamsPerCond = OptionsStruct.ParamsPerCond;
    OptionsStruct = rmfield(OptionsStruct,'ParamsPerCond');
else
    P.ParamsPerCond = {1:P.nParams2Fit};                                    %default: all Params2Fit belong to the first condition
end
P.nConditions = numel(P.ParamsPerCond);

%Set defaults for parameters
DefaultsStruct = BCIsetDefaults(P.nParams2Fit,responsesAVC);
defaultFields = fieldnames(DefaultsStruct);
for i=1:length(defaultFields)
    P.(defaultFields{i,1}) = DefaultsStruct.(defaultFields{i,1});
end

%Overwrite the parameters that we will fit with NaNs (such that we don't bias the initial gridsearch)    
for i=1:P.nParams2Fit
    P.(P.ParamNames2Fit{1,i}) = NaN;
end

%Check the Bounds structure in OptionsStruct
if isfield(OptionsStruct,'Bounds')
    if ~isstruct(OptionsStruct.Bounds)
        error('OptionsStruct.Bounds must be a structure containing a separate field for each of the parameter bounds arrays of size [1x4]');
    end
end

%Check the parallel structure in OptionsStruct
if isfield(OptionsStruct,'parallel')
    if ~isstruct(OptionsStruct.parallel)
        error('OptionsStruct.parallel must be a structure containing one or more of the following fields: BADS, VBMC and/or MCMC');
    end
end

%Check the maxRounds structure in OptionsStruct
if isfield(OptionsStruct,'nAttempts')
    if ~isstruct(OptionsStruct.nAttempts)
        error('OptionsStruct.nAttempts must be a structure containing one or more of the following fields: BADS, VBMC and/or MCMC');
    end
    if isfield(OptionsStruct.nAttempts,'BADS'); OptionsStruct.nAttempts.BADS = max([1 1],OptionsStruct.nAttempts.BADS); end    %Ensure nAttempts [min,max] are 1 or larger
    if isfield(OptionsStruct.nAttempts,'VBMC'); OptionsStruct.nAttempts.VBMC = max([1 1],OptionsStruct.nAttempts.VBMC); end
    if isfield(OptionsStruct.nAttempts,'MCMC'); OptionsStruct.nAttempts.MCMC = max([1 1],OptionsStruct.nAttempts.MCMC); end
end

%Let OptionsStruct overwrite the default values (defined by the user)
optionFields = fieldnames(OptionsStruct);
for i=1:length(optionFields)
    %All subfields separately
    if strcmp(optionFields{i,1},'Bounds') || strcmp(optionFields{i,1},'maxRounds')
        subFields = fieldnames(OptionsStruct.(optionFields{i,1}));
        for j=1:length(subFields)
            P.(optionFields{i,1}).(subFields{j,1}) = OptionsStruct.(optionFields{i,1}).(subFields{j,1});
        end
    %Otherwise the main fields directly    
    else
        P.(optionFields{i,1}) = OptionsStruct.(optionFields{i,1});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check the fixed parameter values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note that fitted parameters will normally be set to NaN, unless the user has provided an initial guess    

if ~isnan(P.Pcommon); P.Pcommon = min(max(0,P.Pcommon),1); end                  % 0 <= Pcommon <= 1
if ~isnan(P.CSJthresh); P.CSJthresh = min(max(0,P.CSJthresh),1); end            % 0 <= CSJthresh <= 1
if ~isnan(P.muP); P.muP = min(max(-90,P.muP),90); end                           % -90 <= muP <= 90

if ~isnan(P.sigmaA); P.sigmaA = min(max(1e-3,P.sigmaA),1000); end               % 1e-3 <= sigmaA <= 1000
if ~isnan(P.sigmaV); P.sigmaV = min(max(1e-3,P.sigmaV),1000); end               % 1e-3 <= sigmaV <= 1000
if ~isnan(P.sigmaP); P.sigmaP = min(max(1,P.sigmaP),1000); end                  % 1 <= sigmaP <= 1000    --> minimum is 1 to avoid numerical errors (see also below for lower bound)
if ~isnan(P.sigmaMotor); P.sigmaMotor = min(max(1e-3,P.sigmaMotor),1000); end   % 1e-3 <= sigmaMotor <= 1000

if ~isnan(P.deltaSigmaA); P.deltaSigmaA = max(0,P.deltaSigmaA); end             % 0 <= deltaSigmaA
if ~isnan(P.deltaSigmaV); P.deltaSigmaV = max(0,P.deltaSigmaV); end             % 0 <= deltaSigmaV
if ~isnan(P.deltaXA); P.deltaXA = max(-1+1e-3,P.deltaXA); end                   % -1+1e-3 <= deltaXA
if ~isnan(P.deltaXV); P.deltaXV = max(-1+1e-3,P.deltaXV); end                   % -1+1e-3 <= deltaXV

%Check that the lapse rates are larger than 1e-9 (to avoid infinite likelihood issues: log(prob=0) = -inf) and lower than 1.  
if ~isnan(P.lapseR); P.lapseR = min(max(1e-9,P.lapseR),1-1e-9); end        
if ~isnan(P.lapseRA); P.lapseRA = min(max(1e-9,P.lapseRA),1-1e-9); end  
if ~isnan(P.lapseRV); P.lapseRV = min(max(1e-9,P.lapseRV),1-1e-9); end  
if ~isnan(P.lapseRC); P.lapseRC = min(max(1e-9,P.lapseRC),1-1e-9); end  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert the bounds %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Collect the bounds for the parameters of interest
P.Bounds.LB = nan(1,P.nParams2Fit);
P.Bounds.PLB = nan(1,P.nParams2Fit);
P.Bounds.PUB = nan(1,P.nParams2Fit);
P.Bounds.UB = nan(1,P.nParams2Fit);
for i=1:P.nParams2Fit
    P.Bounds.LB(1,i) = P.Bounds.(P.ParamNames2Fit{1,i})(1,1);
    P.Bounds.PLB(1,i) = P.Bounds.(P.ParamNames2Fit{1,i})(1,2);
    P.Bounds.PUB(1,i) = P.Bounds.(P.ParamNames2Fit{1,i})(1,3);
    P.Bounds.UB(1,i) = P.Bounds.(P.ParamNames2Fit{1,i})(1,4);
    
    %Check the bounds for the parameters that we fit
    %We don't use 0 and 1 as hard bounds for params that we fit in logit-space, (or 0 for params in log-space)   
    %because these amount to -inf and inf, which my prior function and VBMC don't like
    if any(strcmp(P.ParamNames2Fit{1,i},{'Pcommon','CSJthresh','lapseR','lapseRA','lapseRV','lapseRC'}))                                %These params (between 0 and 1) are fitted in logit space
        if P.Bounds.LB(1,i) < 1e-9
            P.Bounds.LB(1,i) = 1e-9;
            P.Bounds.PLB(1,i) = max(P.Bounds.PLB(1,i),P.Bounds.LB(1,i)+1e-3);
        end
        if P.Bounds.UB(1,i) > (1-1e-9)
            P.Bounds.UB(1,i) = (1-1e-9);
            P.Bounds.PUB(1,i) = min(P.Bounds.PUB(1,i),P.Bounds.UB(1,i)-1e-3);
        end   
    elseif any(strcmp(P.ParamNames2Fit{1,i},{'deltaXA','deltaXV'}))                                                                     %Ensure that 'deltaXA' and 'deltaXV' are postive, to avoid midline crossings...        
        if P.Bounds.LB(1,i) < -1+1e-3                                                                                                   %... and to avoid being equal to -1 exactly which would make all xA exactly muP                
            P.Bounds.LB(1,i) = -1+1e-3;
            P.Bounds.PLB(1,i) = max(P.Bounds.PLB(1,i),P.Bounds.LB(1,i)+1e-3);
        end      
    elseif any(strcmp(P.ParamNames2Fit{1,i},{'sigmaP','sigmaA','sigmaV','sigmaMotor','deltaSigmaA','deltaSigmaV'}))                     %These positive params are fitted in log-space    
        if P.Bounds.LB(1,i) < 1e-3
            P.Bounds.LB(1,i) = 1e-3;
            P.Bounds.PLB(1,i) = max(P.Bounds.PLB(1,i),P.Bounds.LB(1,i)+1e-3);
        end
        if P.Bounds.UB(1,i) > 1000                                                                                                      %Use realistic minima / maxima
            P.Bounds.UB(1,i) = 1000;
            P.Bounds.PUB(1,i) = min(P.Bounds.PUB(1,i),P.Bounds.UB(1,i)-1);
        end 
    elseif strcmp(P.ParamNames2Fit{1,i},'muP')
        if P.Bounds.LB(1,i) < -90
            P.Bounds.LB(1,i) = -90;
            P.Bounds.PLB(1,i) = max(P.Bounds.PLB(1,i),P.Bounds.LB(1,i)+1);
        end
        if P.Bounds.UB(1,i) > 90                                                                                                        %Use realistic minima / maxima
            P.Bounds.UB(1,i) = 90;
            P.Bounds.PUB(1,i) = min(P.Bounds.PUB(1,i),P.Bounds.UB(1,i)-1);
        end  
    end
    
    %Avoid numerical errors. Minimum of sigmaP is 1 instead of 0.001, because if sigmaP and sigmaV are both very small, then all pC become NaN (zero divided by zero) in BCIcompLL_bci.m
    if strcmp(P.ParamNames2Fit{1,i},'sigmaP')
        P.Bounds.LB(1,i) = max(P.Bounds.LB(1,i),1);
        P.Bounds.PLB(1,i) = max(P.Bounds.PLB(1,i),P.Bounds.LB(1,i)+1e-2);
    end
    
    %Check that the bounds are in the right order
    if any([~(P.Bounds.LB(1,i) < P.Bounds.PLB(1,i)), ~(P.Bounds.PLB(1,i) < P.Bounds.PUB(1,i)), ~(P.Bounds.PUB(1,i) < P.Bounds.UB(1,i))])
        error('Check the parameter bounds: they should be in the following order: LB < PLB < PUB < UB')
    end
    
    %Check that initializations of the fitted parameter are within/equal to the plausible bounds - otherwise move them to nearest plausible bound   
    if P.(P.ParamNames2Fit{1,i}) < P.Bounds.PLB(1,i)
        P.(P.ParamNames2Fit{1,i}) = P.Bounds.PLB(1,i);
    elseif P.(P.ParamNames2Fit{1,i}) > P.Bounds.PUB(1,i)
        P.(P.ParamNames2Fit{1,i}) = P.Bounds.PUB(1,i);
    end
end

%Convert the bounds (to log or logit space, in which we fit the parameters)
P.Bounds.conv.LB = BCIconvertParams2Fit(P.Bounds.LB,P,'real2conv');
P.Bounds.conv.PLB = BCIconvertParams2Fit(P.Bounds.PLB,P,'real2conv');
P.Bounds.conv.PUB = BCIconvertParams2Fit(P.Bounds.PUB,P,'real2conv');
P.Bounds.conv.UB = BCIconvertParams2Fit(P.Bounds.UB,P,'real2conv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform some additional checks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check that BADS optimization is performed before VBMC / MCMC 
%Note that this function, BCIsetOptions.m, is not run when BCIfitModel.m is run with previous BCIfitResults as input parameter (in which case BADS may already have run previously and can now be skipped).   
if ~P.BADS && (P.VBMC || P.MCMC)
    warning('VBMC or MCMC was requested without BADS optimization, but these functions require BADS optimized parameters. So we will run BADS optimization first.');
    P.BADS = 1;
end

%Check that VBMC is done before MCMC (necessary to sample good initializations for the MCMC walkers and thus avoid a lenghty burn-in period)
if P.MCMC && ~P.VBMC
    warning('MCMC was requested without VBMC. We will perform VBMC first to speed up MCMC.');
    P.VBMC = 1;
end    

end %[EOF]
