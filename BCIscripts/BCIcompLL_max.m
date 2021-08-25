function [LLvectorA,LLvectorV,LLvectorC] = BCIcompLL_max(P,sA,sV,sA_resp,sV_resp,CSJ_resp)
% Generate response frequencies that maximise the log-likelihood for one condition (combination of sA and sV)
% Generally, this means that we generate distributions that are identical to the actual response distributions.  
% And compute loglike for A, V and CSJ responses

% Gather relevant variables from structure P
nGrid = P.nGridNumInt;
RespLocs = P.RespLocs;
RespRange = P.RespRange;
sigmaMotor = P.sigmaMotor;
triangpdfSpeedUp = P.triangpdfSpeedUp;
integrateMethod = P.integrateMethod;

% Determine number of trials
nTrialsA = numel(sA_resp);
nTrialsV = numel(sV_resp);
nTrialsC = numel(CSJ_resp);

% Initialize output vectors
LLvectorA = nan(nTrialsA,1);
LLvectorV = nan(nTrialsV,1);
LLvectorC = nan(nTrialsC,1);

%Discretized responses
if all(~isnan(RespLocs))
    
    %Create bins to use for the histogram (we use 'histcounts' to avoid problems with numerical approximations)
    nRespLocs = numel(RespLocs);
    if nRespLocs == 1
        edges = [-inf inf]; 
    else
        edges = [-inf mean([RespLocs(1:(end-1)); RespLocs(2:end)],1) inf];  
    end
    
    %Auditory responses
    if (nTrialsA > 0)
        likeFunA = histcounts(RespLocs(sA_resp),edges)/nTrialsA;
        LLvectorA(:) = log(likeFunA(sA_resp));                              %Note that the responses were already transformed to bin_idx numbers before
    end
    
    %Visual responses
    if (nTrialsV > 0)
        likeFunV = histcounts(RespLocs(sV_resp),edges)/nTrialsV;
        LLvectorV(:) = log(likeFunV(sV_resp));                              %Note that the responses were already transformed to bin_idx numbers before
    end

% Default integration method for continuous responses     
elseif integrateMethod == 1
    
    %Determine motorNoisePDF function
    if triangpdfSpeedUp
        motorNoisePDFfun = @(x,mu,sigma) bsxfun_triangpdf(x,mu,sigma);      %Triangular motor noise pdf as an approximation to gaussian motor noise pdf (qtrapz over zeros is fast [for all s_resp far away from s_hat])  
    else
        motorNoisePDFfun = @(x,mu,sigma) bsxfun_normpdf(x,mu,sigma);        %Default Gaussian motor noise pdf
    end
    
    %Auditory responses
    if (nTrialsA > 0)
        sA_hat = sA_resp;                                                   %Create a distribution that is identical to the actual responses (1st dim)
        xA_prob = ones(nTrialsA,1)/nTrialsA;                                %Each response has equal probability  
        
        % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.  
        % For each "response bin" (with centres defined by sA_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sA_hat
        % We can then compute the overall likelihood for each response by integrating out the xAs (1st dim) (i.e. the likelihood is a weighted average pdf) 
        likeA = qtrapz(bsxfun(@times,xA_prob,motorNoisePDFfun(sA_resp', sA_hat, sigmaMotor)),1);
        LLvectorA(:) = log(likeA);                                          %Log transform and save
    end
    
    %Visual responses   
    if (nTrialsV > 0)
        sV_hat = sV_resp';                                                  %Create a distribution that is identical to the actual responses (2nd dim)
        xV_prob = ones(1,nTrialsV)/nTrialsV;                                %Each response has equal probability 
        
        % We assume that observers' responses are imprecise due to motor noise. We use "sigmaMotor" as a parameter that determines the SD of that Gaussian motor noise.  
        % For each "response bin" (with centres defined by sV_hat) the likelihood of any observer's response is given by a Gaussian pdf centred on sV_hat
        % We can then compute the overall likelihood for each response by integrating out the xVs (2nd dim) (i.e. the likelihood is a weighted average pdf)  
        likeV = qtrapz(bsxfun(@times,xV_prob,motorNoisePDFfun(sV_resp, sV_hat, sigmaMotor)),2);
        LLvectorV(:) = log(likeV);                                          %Log transform and save
    end

% Alternative integration method for continuous responses    
elseif integrateMethod == 2
    
    %Auditory responses
    if (nTrialsA > 0)
        sA_hat = sort(sA_resp);                                             %Create a sorted distribution that is identical to the actual responses (1st dim)
        xA_prob = ones(nTrialsA,1)/nTrialsA;                                %Each response has equal probability  
        
        % Integrate out xA using a custom-made version of numerical intergration that takes into account the width of each sA_hat bin in the xA direction
        [likeFunA,likeFunA_grid,likeFunA_spacing] = BCIcompLikeFun(sA_hat,xA_prob,nGrid,RespRange,1,sigmaMotor);
        LLvectorA(:) = log(lininterp1(likeFunA_grid,likeFunA,sA_resp,[],likeFunA_spacing));                 %Interpolate likelihood function to find pdf value for each response 
    end
    
    %Visual responses   
    if (nTrialsV > 0)
        sV_hat = sort(sV_resp)';                                            %Create a sorted distribution that is identical to the actual responses (2nd dim)
        xV_prob = ones(1,nTrialsV)/nTrialsV;                                %Each response has equal probability 
        
        % Integrate out xV using a custom-made version of numerical intergration that takes into account the width of each sV_hat bin in the xV direction    
        [likeFunV,likeFunV_grid,likeFunV_spacing] = BCIcompLikeFun(sV_hat,xV_prob,nGrid,RespRange,2,sigmaMotor);
        LLvectorV(:) = log(lininterp1(likeFunV_grid,likeFunV,sV_resp,[],likeFunV_spacing));                 %Interpolate likelihood function to find pdf value for each response 
    end
end    

%Common source judgments
if (nTrialsC > 0)
    if isnan(sA) || isnan(sV)
        %Unisensory - fixed
        likeFunCSJ = [0 1]; %[Common = 0, Indep = 1];
    else
        %Bisensory - copy frequencies of actual responses
        likeFunCSJ = nan(1,2);
        likeFunCSJ(1) = sum(CSJ_resp == 1)/nTrialsC;
        likeFunCSJ(2) = sum((CSJ_resp ~= 1) & ~isnan(CSJ_resp))/nTrialsC;   
    end
    LLvectorC(:) = log(likeFunCSJ(CSJ_resp));                               %Note that the responses were already transformed to bin_idx numbers before
end
    
end %[EOF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function z = qtrapz(y,dim)
%QTRAPZ  Quick trapezoidal numerical integration.
%   Z = QTRAPZ(Y) computes an approximation of the integral of Y via
%   the trapezoidal method (with unit spacing).  To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, QTRAPZ(Y) is the integral of Y. For matrices, QTRAPZ(Y)
%   is a row vector with the integral over each column. For N-D
%   arrays, QTRAPZ(Y) works across the first non-singleton dimension.
%
%   Z = QTRAPZ(Y,DIM) integrates across dimension DIM of Y. The length of X 
%   must be the same as size(Y,DIM).
%
%   QTRAPZ is up to 3-4 times faster than TRAPZ for large arrays.
%
%   See also TRAPZ.

% Luigi Acerbi <luigi.acerbi@nyu.edu>
% Version 1.0. Release date: Jul/20/2015.

% By default integrate along the first non-singleton dimension
if nargin < 2; dim = find(size(y)~=1,1); end    

% Behaves as sum on empty array
if isempty(y); z = sum(y,dim); return; end

% Compute dimensions of input matrix    
if isvector(y); n = 1; else n = ndims(y); end

switch n
    case {1,2}      % 1-D or 2-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:) + y(end,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1) + y(:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end

    case 3      % 3-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:) + y(end,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:) + y(:,end,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1) + y(:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    case 4      % 4-D array
        switch dim
            case 1
                z = sum(y,1) - 0.5*(y(1,:,:,:) + y(end,:,:,:));
            case 2
                z = sum(y,2) - 0.5*(y(:,1,:,:) + y(:,end,:,:));
            case 3
                z = sum(y,3) - 0.5*(y(:,:,1,:) + y(:,:,end,:));
            case 4
                z = sum(y,4) - 0.5*(y(:,:,:,1) + y(:,:,:,end));
            otherwise
                error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
        end                

    otherwise   % 5-D array or more
        for iDim = 1:n; index{iDim} = 1:size(y,iDim); end
        index1 = index;     index1{dim} = 1;
        indexend = index;   indexend{dim} = size(y,dim);
        try
            z = sum(y,dim) - 0.5*(y(index1{:}) + y(indexend{:}));
        catch
            error('qtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');            
        end
end

end %[EoF]

%-------------------------------------------------------------------------

function y = bsxfun_normpdf(x,mu,sigma)
%BSXFUN_NORMPDF Vectorized normal probability density function (pdf).
%   Y = BSXFUN_NORMPDF(X,MU,SIGMA) returns the pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMCDF, NORMPDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, x - mu, sigma).^2), sigma)/sqrt(2*pi);
    elseif isscalar(sigma)
        y = exp(-0.5*(bsxfun(@minus, x, mu)/sigma).^2)/(sigma*sqrt(2*pi));
    else
        y = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma).^2), sigma)/sqrt(2*pi);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

end %[EoF]

%-------------------------------------------------------------------------

function y = bsxfun_triangpdf(x,mu,sigma)
%BSXFUN_TRIANGPDF Vectorized probability density function (pdf) of
%symmetrical triangular distribution as an approximation to the normal pdf. 
%   Y = BSXFUN_TRIANGPDF(X,MU,SIGMA) returns the pdf of the symmetrical
%   triangular distribution with mean MU. Its width is defined relative to 
%   width parameter SIGMA such that the interval MU +/- SIGMA contains ~68%
%   of its mass (similar to a normal distribution with same mu and sigma).
%   Dimensions of X, MU, and SIGMA must either match, or be equal to one. 
%   Computation of the pdf is performed with singleton expansion enabled 
%   via BSXFUN. The size of Y is the size of the input arguments (expanded
%   to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   Adapted from BSXFUN_NORMPDF by Luigi Acerbi.
%
%   This function is about 3 times faster than BSXFUN_NORMPDF

if nargin<3
    error('bmp:bsxfun_triangpdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

widthFactor = 2.29; %--> Create a symmetrical triangular pdf where 68% of its mass falls within the interval that would also contain 68% of the normal distribution's mass (i.e. interval is mu +/- SD)
                         %This also means that 97.8% of the gaussian mass is captured within the interval of the triangular pdf (outside that interval pdf = 0 for triangular, and pdf = small for gaussian)   

%y = ((widthFactor*sigma)-abs(mu-x))/((widthFactor*sigma)^2)
try
    if isscalar(mu)
        y = bsxfun(@rdivide,bsxfun(@minus,(widthFactor*sigma),abs(mu-x)),(widthFactor*sigma).^2);        
    elseif isscalar(sigma)
        y = ((widthFactor*sigma)-abs(bsxfun(@minus,mu,x)))/((widthFactor*sigma)^2);
    else
        y = bsxfun(@rdivide,bsxfun(@minus,(widthFactor*sigma),abs(bsxfun(@minus,mu,x))),(widthFactor*sigma).^2);
    end
catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end

y = max(y,0);       %This is significantly (!) faster than y(y<0) = 0;

end %[EoF]

%-------------------------------------------------------------------------

function Vout = lininterp1(X,V,Xq,extrap,deltaX)
%Simple 1D linear interpolation on a regular ascending grid.
%Inputs X, V and Xq must be vectors. 
%Output Vout is a column vector.

%This function is a simplified version of the lininterp1 function that was
%originally written by Luigi Acerbi: https://github.com/lacerbi/lautils-mat
%While Luigi's version supports matrices as inputs, this simplified
%version is about 4x faster, and nearly 10x faster than Matlab's interp1.

%Determine length of the given grid 
%Note that X is allowed to be given as X = [min max];
nGrid = length(V);

%Ensure Xq is not empty and a column vector
if isempty(Xq)
    Vout = [];
    return
else
    Xq = Xq(:);
end

%Set default extrapolation: NaN
if nargin < 4
    extrap = NaN;
end

%Determine deltaX from X if not already supplied
if nargin < 5
    deltaX = (X(end)-X(1))/(nGrid-1);
end

%Find Xq out of bounds
flag1 = Xq < X(1);
flag2 = Xq > X(end);
flag3 = isnan(Xq);

%Transform Xq as if X was given as 1:length(X)
Xq = (Xq-X(1))/deltaX + 1;                                                                          
                                
%Indices of nearest lower neighbour X
Xqi = floor(Xq);
Xqi(flag1 | flag2 | flag3) = 1;     %temp --> will be overwritten below

%Distance from Xqi
delta = Xq - Xqi; 

%Expand to avoid Xqi==numel(V) errors when asking for Xqi+1 index below
V = [V(:); 0];

%Distance-weighted average (i.e. linear interpolation)   
Vout = (1-delta).*V(Xqi) + delta.*V(Xqi+1);

%Set values outside of X (extrapolation)
if isempty(extrap)                                              
    %Xq outside X get V values of nearest edge.
    Vout(flag1) = V(1);        
    Vout(flag2) = V(nGrid);
else
    %Xq outside X get specified value (default = NaN)
    Vout(flag1 | flag2) = extrap;
end
Vout(flag3) = NaN; 

end %[EoF]
