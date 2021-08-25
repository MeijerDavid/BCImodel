function paramsConverted = BCIconvertParams2Fit(params,P,direction)
%Convert parameters to log/logit space or back (depending on the
%'direction': either 'real2conv' or 'conv2real'). If 'params' is a matrix,
%then conversion is performed on the entire column.

%Initialize
paramsConverted = params;

%Convert
if strcmp(direction,'real2conv')
    for i=1:P.nParams2Fit
        if any(strcmp(P.ParamNames2Fit{1,i},{'Pcommon','CSJthresh','lapseR','lapseRA','lapseRV','lapseRC'}))                                %Fit params between 0 and 1 in logit space
            paramsConverted(:,i) = logit(params(:,i));
        elseif any(strcmp(P.ParamNames2Fit{1,i},{'sigmaP','sigmaA','sigmaV','sigmaMotor','deltaSigmaA','deltaSigmaV'}))                     %Fit positive params in log-space
            paramsConverted(:,i) = log(params(:,i));
        elseif any(strcmp(P.ParamNames2Fit{1,i},{'deltaXA','deltaXV'}))                                                                     %Ensure positivity (LB = 1e-9), then fit in log-space                                            
            paramsConverted(:,i) = log(params(:,i)-P.Bounds.(P.ParamNames2Fit{1,i})(1)+1e-9);    
        end
    end

%Convert back    
elseif strcmp(direction,'conv2real')
    for i=1:P.nParams2Fit
        if any(strcmp(P.ParamNames2Fit{1,i},{'Pcommon','CSJthresh','lapseR','lapseRA','lapseRV','lapseRC'}))                                %Inverse of logit is logistic
            paramsConverted(:,i) = logistic(params(:,i));
        elseif any(strcmp(P.ParamNames2Fit{1,i},{'sigmaP','sigmaA','sigmaV','sigmaMotor','deltaSigmaA','deltaSigmaV'}))                     %Inverse of log is exponential
            paramsConverted(:,i) = exp(params(:,i));
        elseif any(strcmp(P.ParamNames2Fit{1,i},{'deltaXA','deltaXV'}))                                                                     %Inverse of log is exponential, then return the shift 
            paramsConverted(:,i) = exp(params(:,i))+P.Bounds.(P.ParamNames2Fit{1,i})(1)-1e-9;    
        end
    end

%Error    
else 
    error('Direction of parameter conversion is UNKNOWN');
end

end %[EOF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function p = logistic(alpha)

p = 1./(1 + exp(-alpha));

end %[EoF]

% ------------------------------------------------------------------------

function alpha = logit(p)

alpha = log(p./(1-p));

end %[EoF]
