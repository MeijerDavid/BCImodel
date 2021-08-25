function defaults = BCIsetDefaults(nParams2Fit,responsesAVC)
% Define the defaults parameter settings and parameter bounds. 
% This function is called by "BCIsetOptions.m"

%Parameters that can be fitted
defaults.Pcommon = 0.5;                     %common source prior
defaults.sigmaA = 5;                        %sigma of auditory noise
defaults.sigmaV = 2;                        %sigma of visual noise
defaults.sigmaP = 10;                       %sigma of spatial prior
defaults.muP = 0;                           %centre of spatial prior

defaults.deltaSigmaA = 0;                   %Linear sigmaA shift with eccentricity: sigmaA_final = sigmaA + deltaSigmaA*abs(sA-muP)  
defaults.deltaSigmaV = 0;                   %Linear sigmaV shift with eccentricity: sigmaV_final = sigmaV + deltaSigmaV*abs(sV-muP)  
defaults.deltaXA = 0;                       %Linear xA shift with eccentricity: xA_final = xA + deltaXA*xA 
defaults.deltaXV = 0;                       %Linear xV shift with eccentricity: xV_final = xV + deltaXV*xV --> this creates non-Gaussian likelihood functions: see Wei & Stocker (2015) Nature Neuroscience

defaults.CSJthresh = 0.5;                   %common source judgment threshold 
defaults.sigmaMotor = 1;                    %sigma of motor noise for continuous responses

defaults.lapseR = NaN;                      %lapse rate for all responses (if ~isnan(lapseR), then lapseR overwrites lapseRA, lapseRV and lapseRC; see BCIcomputeLL.m)
defaults.lapseRA = 0;                       %lapse rate for auditory responses only
defaults.lapseRV = 0;                       %lapse rate for visual responses only   
defaults.lapseRC = 0;                       %lapse rate for common-source judgments

%Other settings
defaults.ForceLogLikelihood = 0;            %Should we use log-likelihood instead of log-posterior probability?     
defaults.nGridSearch = 1000*nParams2Fit;    %Nr of gridpoints for gridsearch that precedes BADS
defaults.nGridNumInt = 101;                 %Nr of gridpoints in each dimension that is used for numerical integration over all internal representations xA or xV
defaults.decisionFun = 'ModelAveraging';    %BCI decision function (see Wozny, Beierholm, Shams, 2010, PLOS Comp. Biology)

defaults.RespLocs = NaN;                    %Set vector of response locations for button responses. Leave as NaN for continuous responses (e.g. mouse responses) 

%Determine possible response range for continuous responses based on the given responses    
defaults.RespRange = [min(responsesAVC(:,[1 2]),[],'all'), max(responsesAVC(:,[1 2]),[],'all')];

%Additional options for continuous responses   
defaults.integrateMethod = 1;               %Integrate method 2 is slightly less precise but faster when there are many responses per location (computation time increases linearly with number of responses per loc in method 1!)
defaults.triangpdfSpeedUp = 0;              %Using a triangular pdf approximation instead of the true normpdf is about 3 times faster!

%Set default parallel computing behavior (logical switches)
defaults.parallel.BADS = 0;                 
defaults.parallel.VBMC = 0;                 
defaults.parallel.MCMC = 0;                 

%Set default number of attempts [min,max] to reach convergence and select the best result    
defaults.nAttempts.BADS = [1 4];              
defaults.nAttempts.VBMC = [1 4];             
defaults.nAttempts.MCMC = [1 4];                        

%Computational algorithms to run
defaults.BADS = 1;                          %Switch for BADS. Note: some other BADS default settings may be set in 'BCIbads.m'
defaults.VBMC = 0;                          %Switch for VBMC. Note: some other VBMC default settings may be set in 'BCIvbmc.m'
defaults.MCMC = 0;                          %Switch for MCMC. Note: some other MCMC default settings may be set in 'BCImcmc.m'

%Display results?
defaults.Display = 0;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter bounds %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Please note that the parameter bounds are extremely important as they 
% also define the prior probability distributions! (The priors will be 
% defined as trapezoids with highest probability between the plausible 
% bounds, and linearly decreasing on either side towards the hard bounds). 

%Define default bounds:       [Hard Lower, Plausible Lower, Plausible Upper, Hard Upper]    
defaults.Bounds.Pcommon     = [   0            0.01             0.99             1     ];           %Be realistic when changing the parameter bounds
defaults.Bounds.sigmaA      = [   0            0.1             15              100     ];           %Extraordinarily small/large bounds do not help the fitting algorithms! 
defaults.Bounds.sigmaV      = [   0            0.1             15              100     ];        
defaults.Bounds.sigmaP      = [   0            5              500             1000     ];
defaults.Bounds.muP         = [ -45           -5                5               45     ];

defaults.Bounds.deltaSigmaA = [   0            1e-6             0.5             10     ];           
defaults.Bounds.deltaSigmaV = [   0            1e-6             0.5             10     ];  
defaults.Bounds.deltaXA     = [  -1           -0.5              0.5             10     ];           %Setting a lower bound of -1 ensures that likelihood distributions can move towards the midline, but cannot cross it
defaults.Bounds.deltaXV     = [  -1           -0.5              0.5             10     ]; 

defaults.Bounds.CSJthresh   = [   0            0.1              0.9              1     ];
defaults.Bounds.sigmaMotor  = [   0            0.1             15              100     ];

defaults.Bounds.lapseR      = [   0            1e-6             0.25             1     ];        
defaults.Bounds.lapseRA     = [   0            1e-6             0.25             1     ];        
defaults.Bounds.lapseRV     = [   0            1e-6             0.25             1     ];
defaults.Bounds.lapseRC     = [   0            1e-6             0.25             1     ];

%NOTE: Although we set the default bounds here, some are changed slightly in "BCIsetOptions.m"   
%Those changes are to avoid numerical problems. E.g. We use '1e-9' as hard lower bound instead of '0'.   

end %[EOF]
