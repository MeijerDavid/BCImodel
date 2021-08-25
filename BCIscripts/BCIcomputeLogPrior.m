function lpp = BCIcomputeLogPrior(params,P)
%Compute the log prior probability of the parameters
%We construct the prior PDFs as trapezoids (piecewise continuous in the converted domain): 
%They are maximal between PLB and PUB, and linearly going to zero on either side: from PLB to LB, and from PUB to UB. 

%Gather bounds from struct P (note that these bounds are set in BCIsetParams.m and are not changed afterwards)    
LB = P.Bounds.conv.LB;
UB = P.Bounds.conv.UB;
PLB = P.Bounds.conv.PLB;
PUB = P.Bounds.conv.PUB;

%Number of parameters
nParams = numel(params);

%Compute the height such that the area under the prior integrates to 1.
maxHeightPriors = 1./(0.5*(PLB-LB)+(PUB-PLB)+0.5*(UB-PUB));

%Compute the prior probability per parameter
pp = nan(1,nParams);
for i=1:nParams
    
    %Outside/on hard bounds
    if (params(i) <= LB(i)) || (params(i) >= UB(i))
        pp(i) = 0;
    %Within plausible bounds    
    elseif (params(i) >= PLB(i)) && (params(i) <= PUB(i))
        pp(i) = maxHeightPriors(i);
    %On upward slope
    elseif (params(i) > LB(i)) && (params(i) < PLB(i))
        extend = (params(i)-LB(i))/(PLB(i)-LB(i));
        pp(i) = extend*maxHeightPriors(i);
    %On downward slope 
    elseif (params(i) < UB(i)) && (params(i) > PUB(i))
        extend = (params(i)-PUB(i))/(UB(i)-PUB(i));
        pp(i) = (1-extend)*maxHeightPriors(i);
    end
end

%Avoid log(0) = -inf issues
pp(pp < realmin) = realmin;

%Compute the overall log prior probability
lpp = sum(log(pp));

end %[EOF]