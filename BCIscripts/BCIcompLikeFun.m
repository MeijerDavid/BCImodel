function [likeFun,likeFun_grid,likeFun_spacing] = BCIcompLikeFun(s_hat,prob,nGrid,RespRange,interpDim,sigmaMotor)

%Ensure a relevant minimum grid-spacing (relative to sigmaMotor) --> for reasons see below   
min_spacing = sigmaMotor/100;

%Prepare a regular grid (with spacing > min_spacing) that captures all of the s_hat samples  
min_grid = max(min(s_hat(:)),RespRange(1));
max_grid = min(max(s_hat(:)),RespRange(2));
if min_grid < max_grid %Normally
    suggested_spacing = (max_grid-min_grid)/(nGrid-1);
    if suggested_spacing < min_spacing                                      %If suggested spacing is too small, then choose fewer grid-points on the same interval such that the grid-spacing has practical relevance
        nGrid = max(2,floor((max_grid-min_grid)/min_spacing));
    end
else %min_grid >= max_grid --> Uncommon case                                                     
    spacing2coverFullRange = (RespRange(2)-RespRange(1))/(nGrid-1);
    if max_grid <= RespRange(1)                                             %If all s_hat <= RespRange(1)
        max_grid = min_grid + spacing2coverFullRange;
    elseif min_grid >= RespRange(2)                                         %If all s_hat >= RespRange(2)    
        min_grid = max_grid - spacing2coverFullRange;
    end
    nGrid = 2; %Perform integration on these 2 gridpoints only, but the grid will be padded to the full response range before the Gaussian convolution below..
end
likeFun_grid = linspace(min_grid,max_grid,nGrid);
likeFun_spacing = (likeFun_grid(nGrid)-likeFun_grid(1))/(nGrid-1);

%Integrate the s_hat and accompanying probabilities on the regular grid
precisionFactor = 10;                                                                                                   %Default precision factor for the Riemann Sum
likeFun = numIntegrAndInterp(s_hat,prob,likeFun_grid,interpDim,precisionFactor,likeFun_spacing)/likeFun_spacing;        %Also transform probabilities into probability densities

%Zero-pad up to response range (allow for a wider grid to accomodate motor noise, see below)
PadBegin = fliplr((likeFun_grid(1)-likeFun_spacing):-likeFun_spacing:RespRange(1));
PadEnd = (likeFun_grid(nGrid)+likeFun_spacing):likeFun_spacing:RespRange(2);            %If the grid-spacing is too small, then the padded grid may have (extremely) many grid-points
likeFun_grid = [PadBegin likeFun_grid PadEnd];                                          %Moreover, the length of the motor-noise kernel depends on the grid-spacing too, so the convolution may take extremely long to compute   
likeFun = [zeros(1,numel(PadBegin)) likeFun zeros(1,numel(PadEnd))];                    %That is why we set the minimum (relevant) spacing above

%Convolve with Gaussian kernel to implement motor noise
likeFun = convGauss(likeFun,likeFun_spacing,sigmaMotor);                             

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function pVec = numIntegrAndInterp(xMat,pMat,xVec,interpDim,precisionFactor,xVecRegSpacing)
%Numerically integrate probability bins pMat on accompanying irregular and 
%unsorted grid xMat, and interpolate to form a probability distribution on
%the strictly ascending vector xVec. Both xMat and pMat may be matrices: 
%the operational dimension over which integration and interpolation take 
%place can be specified with interpDim (default = 1). 
%"numIntegrAndInterp" integrates by means of a Riemann Sum over xVec: The 
%probability (pMat) is divided across all xVec bins that are covered by the
%original xMat bins. The precision of this probability division can be set
%by precisionFactor (default = 10; scalar: increase for greater precision). 
%A small speed-up can be achieved if xVec is regularly spaced. If so,
%specify the xVec bin-width with xVecRegSpacing (scalar; default = false,
%meaning that xVec is assumed to have irregular intervals). 

%Set Defaults
if nargin < 4
    interpDim = 1;              %Direction of binning and interpolation on xMat and pMat
end                            
if nargin < 5
    precisionFactor = 10;       %Larger factors increase precision of interpolation but take slightly longer to compute and require more memory. Ensure that precisionFactor is a positive integer.
end                             
if nargin < 6                   
    xVecRegSpacing = false;     %Set to non-false if xVec is regularly spaced, where xVecRegSpacing is the regular 'bin-width'. Providing this information leads to a small speed-up but is not necessary.
end                             %Default assumes that xVec is not regularly spaced ('0')

%Check that the dimension over which to integrate is non-singleton
if size(xMat,interpDim) == 1
    pVec = 0*xVec;              %Return an integral of zeros (similar to 'trapz')
    return
end

%Interpolation and integration is performed over the first dimension of xMat
%If it is requested over another dimension, then make that dimension the first dimension   
nDims = numel(size(xMat));
if interpDim~=1
    if nDims <= 2
        xMat = permute(xMat,[2 1]);                     %Simple transpose
        pMat = permute(pMat,[2 1]);
    elseif nDims > 2
        OtherDims = setdiff(1:nDims,interpDim);         
        xMat = permute(xMat,[interpDim OtherDims]);
        pMat = permute(pMat,[interpDim OtherDims]);
    end
end

%Ensure 2D matrices for xMat and pMat (reshape if necessary - i.e. concatenate all dimensions > 2 into the second dimension)
[nRows,nCols] = size(xMat);
if nDims > 2
    xMat = reshape(xMat,[nRows nCols]);
    pMat = reshape(pMat,[nRows nCols]);
end

%Ensure column vector for xVec
if isrow(xVec)
    xVec = xVec';       
    xVecRowBool = 1;
else
    xVecRowBool = 0;
end
nGrid = length(xVec);

%Create a reference vector for the small xVec bins to the larger original xVec bins
nBins = nGrid*precisionFactor;
if precisionFactor ~= 1
    small2big = zeros(nBins,1);
    small2big(1:precisionFactor:end) = 1;
    small2big = cumsum(small2big);                                          %To be used in call to accumarray further below
end

%Discretize xMat
if xVecRegSpacing
    xVecRegSpacing = xVecRegSpacing/precisionFactor;
    iBins = min(nBins,max(1,round((xMat-xVec(1))/xVecRegSpacing)+1));       %Easy and fast, making use of the regular grid with binWidths defined by xSpacing
else                                                                        
    %Use Matlab's 'discretize' function
    if nGrid == 1
        xVec_edges = [xVec(1)-1; xVec(1)+1];                                %If xVec only has one value, then assume arbitrary stepsize of 1 (choice does not matter)
    else    
        xVec_edges = [xVec(1)-diff(xVec(1:2))/2; xVec(1:(nGrid-1))+diff(xVec)/2; xVec(nGrid)+diff(xVec((nGrid-1):nGrid))/2];
    end
    if precisionFactor ~= 1
        xVec_edges = lininterp1(1:(nGrid+1),xVec_edges,linspace(1,nGrid+1,nBins+1)',NaN,1);   %Break each bin into precisionFactor number of equally sized smaller bins     
    end
    xVec_edges([1 end]) = [-inf inf];                                       %Ensure that the bins at the edges cover all xMat values
    iBins = discretize(xMat,xVec_edges);
end

%How many bins are covered with each xMat step?
bins_cov = reshape(abs(diff(iBins))+1,[(nRows-1)*nCols 1]);

%Find minimum and maximum bin indices for each xMat step (prepare for a simple cumsum further below)
iBin_min = reshape(min(iBins(1:(nRows-1),:),iBins(2:nRows,:)),[(nRows-1)*nCols 1]);
iBin_max = iBin_min + bins_cov - 1;

%Prepare to distribute the pMat probabilities across each of the xVec bins
pMat = filter([0.5 0.5],1,pMat);                                            %Moving average filter of width 2: i.e. pNew(i) = 0.5*(pOld(i-1)+pOld(i))
pBins = reshape(pMat(2:nRows,:),[(nRows-1)*nCols 1])./bins_cov;             %Note: half the probs at each edge went 'missing' (similar to 'trapz')
    
%The suggestion for the following crucial line came from: https://uk.mathworks.com/matlabcentral/answers/69413-how-many-times-does-each-number-appear (Teja Muppirala)      
pVec = cumsum(accumarray(iBin_min,pBins,[nBins+1 1]) + accumarray(iBin_max+1,-pBins,[nBins+1 1]));
pVec(nBins+1) = []; 

%Sum probabilities of the small bins that belong to one large bin
if precisionFactor ~= 1
    pVec = accumarray(small2big,pVec,[nGrid 1]); 
end

%Ensure output vector has same dimension as input vector
if xVecRowBool
    pVec = pVec';
end

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

%-------------------------------------------------------------------------

function y_conv = convGauss(y,deltaX,sigma,alpha)
%convolve a signal "y" with a gaussian kernel. Signal y must be based on a
%regular grid with stepsize "deltaX". The Gaussian kernel will have an SD
%equal to "sigma". The width/order of the kernel is further defined by 
%"alpha", which is the maximally allowed value in terms of +/-SD. 

%Set defaults
if nargin < 2
    deltaX = 1;    
end
if nargin < 3
    sigma = 1;    
end
if nargin < 4
    alpha = 2.5;    %equal to the default of "gausswin"
end

%Ensure row vector y
if size(y,1) > 1
    y = y';
end
nY = size(y,2);

%Create a Gaussian kernel with the correct width 
xTmp = deltaX:deltaX:(alpha*sigma);                                         %Kernel order is defined by alpha*sigma relative to deltaX
xTmp = [fliplr(xTmp) 0 xTmp];                                               %Minimum kernel order is 1 (if deltaX > alpha*sigma) and it's always an odd integer
kernelOrder = numel(xTmp);
gaussKernel = normpdf(xTmp,0,sigma);
gaussKernel = gaussKernel/sum(gaussKernel);                                 %Normalize pdf to probabilities

%Pad the input signal 
nPad = (kernelOrder-1)/2;                                                   
if nY > nPad
    y_pad = [y((1+nPad):-1:2), y, y((nY-1):-1:(nY-nPad))];                  %Reflection padding (least disturbance to signal) --> default
elseif nY > 1
    nPad1 = nY-1;
    nPad2 = nPad-nPad1;
    y_pad = [zeros(1,nPad2), y((1+nPad1):-1:2), y, y((nY-1):-1:(nY-nPad1)),zeros(1,nPad2)]; %Partial reflection padding and zeros padding
else
    y_pad = [y*ones(1,nPad), y, y*ones(1,nPad)];                            %Replication padding if signal consists of only one value (not a very sensible convolution..)
end
%y_pad = [y(1)*ones(1,nPad), y, y(end)*ones(1,nPad)];                       %Replication padding (gives greater weight to values at the edges)  
%y_pad = [zeros(1,nPad), y, zeros(1,nPad)];                                 %Zero padding (gives less weight to values at the edges)  

%Perform the convolution using the standard Matlab function
y_conv = conv(y_pad,gaussKernel,'valid');

end %[EoF]
