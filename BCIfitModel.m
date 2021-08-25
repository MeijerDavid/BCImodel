function BCIfitResults = BCIfitModel(trueLocsAV,responsesAVC,i_Conditions,OptionsStruct)
% BCIfitResults = BCIfitModel(trueLocsAV,responsesAV,i_Conditions,OptionsStruct)
%
% Fit the Bayesian Causal Inference model (Körding, Beierholm, Ma, Quartz, 
% Tenenbaum, Shams, 2007, PLOS One) to experimental data of one subject.  
%
% Alternatively call as: BCIfitResults_new = BCIfitModel(BCIfitResults_old)
% This might be convenient if choosing to run VBMC, MCMC or Display at a 
% later time.
%
% For example usage, please see 'BCIparameterRecoveryScript.m'
%
% INPUTS
%
% *trueLocsAV
% A matrix of size [nTrials x 2]. The first column contains the true 
% locations of the auditory stimuli. The second column contains the true
% locations of the visual stimuli. Use NaN for the missing sensory modality
% in case of unisensory stimuli. 
%
% *responsesAV
% A matrix of size [nTrials x 2] or [nTrials x 3]. The first column 
% contains the responded locations of the auditory stimuli. The second 
% column contains the responded locations of the visual stimuli. The third
% column (optional) contains common source responses: "1" or "2" sources. 
% Use NaN for missing responses. The order of trials should correspond to
% "trueLocsAV".
%
% *i_Conditions (optional)
% A matrix of size [nTrials x nConditions], where each column corresponds
% to an experimental condition. Each column (e.g. i'th) should contain a
% vector of logicals where a 'true' at the j'th row indicates that the j'th
% trial belongs to the i'th condition. One trial may belong to multiple 
% conditions. The order of trials should correspond to "trueLocsAV". See 
% also "ParamsPerCond" further below. Default = true(nTrials,1); i.e. all
% trials belong to the first (and only) condition.
%
% *OptionsStruct (optional)
% Struct with various fields that are used as options in the BCI model:
%
% - "ParamNames2Fit"
%   A cell array of size [1 x nParams2Fit] with the names of the parameters 
%   to fit (one per cell). Default = {'Pcommon','sigmaP','sigmaA','sigmaV'}
%   
%   Choose from the following:  (Default ; [ LB    PLB     PUB     UB ])
%   - "Pcommon"                 ( 0.5    ; [  0    0.01    0.99     1 ])
%   - "sigmaA"                  ( 5      ; [  0    0.1    15      100 ])
%   - "sigmaV"                  ( 2      ; [  0    0.1    15      100 ])
%   - "sigmaP"                  (10      ; [  0    5     500     1000 ])
%   - "muP"                     ( 0      ; [-45   -5       5       45 ])
%   - "lapseR"                  ( 0      ; [  0    1e-6    0.25     1 ]) 
%   - "CSJthresh"               ( 0.5    ; [  0    0.1     0.9      1 ])
%   - "sigmaMotor"              ( 1      ; [  0    0.1    15      100 ])
%   - "deltaSigmaA"             ( 0      ; [  0    1e-6    0.5     10 ]) 
%   - "deltaSigmaV"             ( 0      ; [  0    1e-6    0.5     10 ])  
%   - "deltaXA"                 ( 0      ; [  -1  -0.5     0.5     10 ]) 
%   - "deltaXV"                 ( 0      ; [  -1  -0.5     0.5     10 ])
%
%   where 'muP' is the centre of the spatial prior, 'lapseR' is the lapse
%   rate for each response, 'CSJthresh' is the common source threshold for
%   the 'ModelAveraging' decision function (C=1 if p(C=1) >= CSJthresh), 
%   and 'sigmaMotor' is the motor noise for continuous responses (e.g. 
%   mouse responses). 'deltaSigma*' is an eccentricity-dependent multiplier
%   for the sensory noise: e.g. sigmaA = sigmaA + deltaSigmaA*abs(sA). 
%   Finally, 'deltaX*' is an eccentricity-dependent likelihood shift for
%   the internal spatial representations x: e.g. xA = xA+deltaXA*(xA-muP). 
%   See also: Odegaard, Wozny & Shams (2015) PLoS Comput. Biol.
%
%   Any of the above can instead be fixed parameters (they will not be fit 
%   if unmentioned in "ParamNames2Fit"). In case of fixed parameters they
%   will default to the values in parentheses. If you wish to assign
%   different fixed values to them, then add them as separate fields to 
%   OptionsStruct (e.g. OptionsStruct.Pcommon = 1 for Forced Fusion). In 
%   case the default value for a fitted parameter (i.e. mentioned in 
%   ParamNames2Fit) is changed through OptionsStruct, then the assigned
%   value serves as a best guess, and the grid search (preceding BADS 
%   optimization) will be jittered around this best guess. If no best guess
%   is supplied by the user, then the grid will be sampled from a uniform
%   distribution between the plausible bounds (see below).   
%
%   The parameter bounds (hard lower = LB, plausible lower = PLB,
%   plausible upper = PUB, hard upper = UB) that will be used for fitting 
%   are also given in the parentheses. Any of these can be changed by
%   adding them to OptionsStruct.Bounds as vectors of size [1 x 4]. 
%   E.g. OptionsStruct.Bounds.Pcommon = [0.1 0.3 0.7 0.9]. Beware that the
%   following order should be preserved: [LB < PLB < PUB < UB]. Note
%   further that the choice of bounds is very important because they are
%   also used to define prior probability distributions for the parameters.
%   The priors will be defined as trapezoids with highest probability 
%   between the plausible bounds, and linearly decreasing on either side 
%   towards the hard bounds. See further below under "ForceLogLikelihood".
%
%   If you wish to fit separate lapse rates for the auditory, visual 
%   and/or common-source responses, then add "lapseRA", "lapseRV" and/or 
%   "lapseRC" to ParamNames2Fit (default values and bounds can be changed
%   as above). Alternatively, you can fix these lapse rates by adding them
%   to OptionStruct as: e.g. OptionStruct.lapseRA = 0.05. 
% 
% - "ParamsPerCond"
%   A cell array of size [1 x nConditions]. Each cell should contain a
%   vector of size [1 x nParams2FitInThisCondition]. The vector in cell i
%   contains the indices of the ParamNames2Fit that belong to the i'th
%   experimental condition. E.g. ParamsPerCond = {[1 2 3],[1 2 4]}. 
%   Default = {1:nParams2Fit}, i.e. all params are for the first condition
%   (any other conditions in 'i_Conditions' will be ignored!).
%
% - "ForceLogLikelihood"
%   A boolean (0 = 'no', 1 = 'yes'). By default the BCI model is fit by 
%   maximising log posterior probabilities (log-likelihood + log-prior). 
%   In other words, the resulting parameters correspond to the maximum 
%   a-posteriori (MAP) estimate. If you instead wish to fit the model using
%   maximum likelihood estimation (MLE), then set ForceLogLikelihood = 1. 
%   This option does not affect VBMC/MCMC, which are always fitted using 
%   log-posterior probabilities. Default = 0.
%
% - "RespLocs"
%   A row vector of size [1 x N] where N is the number of different
%   response buttons. Set RespLocs = NaN in case of continuous (e.g. mouse)
%   responses. Default = NaN.
%
% - "RespRange"
%   A row vector of size [1 x 2] that indicates the minimum and maximum,
%   respectively, of the response range (e.g. outer bounds of screen for 
%   mouse responses). This range is used for the uniform distribution from
%   which random responses are sampled (e.g. lapses). This setting is only
%   relevant for continuous responses and not used for discrete responses. 
%   Default = [min(location responses A&V), max(location responses A&V)].
%
% - "decisionFun"
%   String with decision function. Choose either: 'ModelAveraging', 
%   'ModelSelection', or 'ProbabilityMatching'. See Wozny, Beierholm & 
%   Shams (2010) PLoS Comput. Biol. Default = 'ModelAveraging'. 
%
% - "nGridNumInt"
%   A scalar. Number of points on the grid for numerical integration over
%   all internal representations xA or xV. For bisensory trials the total
%   number of samples is equal to nGridNumInt^2. The grid will be centred 
%   on sA,sV and covers the interval between -3 and +3 SD. Default = 101.
%
% - "integrateMethod"
%   A scalar: 1 or 2 (default = 1). This defines the numerical integration
%   method that is used for continuous responses (it does not affect
%   discrete responses). Method 1 computes the likelihood for each response
%   as a probability-weighted average pdf value from the motor noise 
%   distributions centred on the predicted responses. While precise, 
%   computations with method 1 can be rather slow when the data contains 
%   many responses per location combination (sA,sV). Instead, method 2 
%   computes an estimate of the full response distribution (per sA,sV), 
%   convolves it with the motor noise pdf, and interpolates the pdf at the
%   values of the actual responses to obtain the likelihoods. 
%   Alternatively, method 1 can be speeded up (by a factor 2-3) by using a
%   triangular approximation of the Gaussian motor noise pdf. If desired, 
%   set "OptionsStruct.triangpdfSpeedUp = 1" (default = 0).
% 
% - "nGridSearch"
%   A scalar. Before fitting parameters using the optimization function
%   (BADS) we perform a grid search. The grid is formed by choosing at 
%   random nGridSearch values for each param2fit from a uniform 
%   distribution that is bounded between the plausible bounds. 
%   Alternatively, if a best guess is provided for a fitted parameter
%   (e.g. OptionsStruct.Pcommon = 0.5 while OptionsStruct.ParamNames2Fit
%   contains "Pcommon") then the grid is sampled from a Gaussian that is 
%   centred on the best guess and with SD = (PUB-PLB)/6. 
%   Default = 1000*nParams2Fit;
%
% - "nAttempts"
%   A structure with three fields: 'BADS','VBMC' and 'MCMC' (all optional).
%   Each entry must be a [1 x 2] vector indicating the minimum and maximum
%   number of attempts that the algorithm will run to try to converge. 
%   Default = [1 4] for each. For BADS and VBMC this means that each
%   algorithm is run at least 'min' times and at most 'max' times with 
%   different initializations on each run. The best solution is returned.
%   For MCMC, any new round means a continuation of the previous round with
%   an increased thinning level: thinning factor = 2^(roundX-1). I.e. for
%   the default nAttempts.MCMC = [1 4], the maximum thinning factor = 8.     
%
% - "parallel" (memory intensive!) 
%   A structure with three fields: 'BADS','VBMC' and 'MCMC' (all optional).
%   Either of the three logical switches can be set to 1 (default = 0) if
%   one wants to make use of parallel computing (using functions from the 
%   respective Matlab Toolbox). BADS and VBMC will run nParallelWorkers 
%   attempts in parallel on the available CPU cores. While available, 
%   setting "OptionsStruct.parallel.MCMC = 1" is not recommended because 
%   the implemented method is invalid: Each CPU core works on all MCMC 
%   chains simultaneously (differently jittered initializations per core) 
%   and their results are then concatenated to form the full chains. This 
%   implementation breaks the Markov property of the chains. Instead, call
%   "BCIfitResults = BCIfitModel(BCIfitResults)" in parallel for multiple 
%   subjects with BCIfitResults.settings.BADS = 0, VBMC = 0 and MCMC = 1 
%   (e.g. see "BCImcmcGroupPar").  
%
% - "BADS" (default fitting algorithm)
%   Boolean for whether we need to perform Bayesian Adaptive Direct Search
%   (1='yes' or 0='no'). This is the default fitting algorithm
%   for BCIfitModel. BADS is a parameter optimization algorithm by Luigi 
%   Acerbi and Wei Ji Ma (https://github.com/lacerbi/bads/; Acerbi & Ma, 
%   2017, Advances in NIPS). Default = 1. Do not change this option on any
%   first call to "BCIfitModel". If, after optimization of the parameters 
%   you call 'BCIfitModel.m' again in order to subsequently perform VBMC or 
%   MCMC (use 'BCIfitResults' as the only input argument) you should set 
%   BCIfitResults.settings.BADS = 0 in order to avoid re-fitting the
%   parameters using BADS. 
%
% - "VBMC" (rough but fast approximation of full posterior)
%   Boolean for whether we need to perform Variational Bayesian Monte Carlo
%   (1='yes' or 0='no'). This option makes use of the 'vbmc' function by
%   Luigi Acerbi: https://github.com/lacerbi/vbmc. See: Acerbi, 2018, 
%   Advances in NIPS. Default = 0, unless MCMC has been requested, see 
%   below. In that case VBMC is run by default to obtain an approximation 
%   of the posterior that is used to guide/speed-up the MCMC algorithm. 
%   VBMC can be performed at a later time by calling 'BCIfitModel.m'
%   with the output structure of this function ('BCIfitResults') as the 
%   only input argument. In that case, avoid refitting BADS by setting
%   BCIfitResults.settings.BADS = 0.
%
% - "MCMC" (proper but slow approximation of full posterior)
%   Boolean for whether we need to perform Markov Chain Monte-Carlo (MCMC)
%   (1='yes' or 0='no'). Default = 0. The MCMC algorithm implements 
%   parallel slice sampling by making use of the 'eissample_lite' function
%   by Luigi Acerbi (https://github.com/lacerbi/eissample) as was also used
%   in Acerbi, Dokka, Angelaki & Ma, 2018, PLOS Comp. Biology. (N.B. This
%   algorithm does not require a Burn-in period if initializations are 
%   chosen properly: e.g. here they are based on the VBMC approximation). 
%   Thinning level increases adaptively until convergence is reached or 
%   maxRounds.MCMC is exceeded (thinning level T at round X is T = 2^(X-1); 
%   meaning that every Tth sample is saved, all others are discarded). 
%   MCMC can also be performed at a later time by calling 'BCIfitModel.m'
%   with the output structure of this function ('BCIfitResults') as the 
%   only input argument. In that case, avoid refitting BADS/VBMC by setting
%   BCIfitResults.settings.BADS/VBMC = 0. A previous session of MCMC can 
%   also be continued with more thinning by increasing the setting of 
%   'BCIfitResults.settings.maxRounds.MCMC' and calling BCIfitModel again.
%
%   If MCMC is requested, then the following option can be set (default): 
%   - "MCMCcount"       (10000)
%       This is the number of MCMC samples that will be collected. The 
%       actual number of function evaluations (i.e. log-likelihood/prior
%       computations) is likely much higher because the MCMC algorithm does
%       not accept all proposed samples and may apply thinning. E.g. set 
%       OptionsStruct.MCMCcount = 20000 for more MCMC samples but increased 
%       computation time and memory usage.
%
% - "Display" 
%   A boolean (0 = 'no', 1 = 'yes'). If set to 'yes' the program produces
%   some figures that enable comparison of the true response distributions
%   with the BCI-predicted response distributions (after BADS parameter 
%   optimization). For VBMC and MCMC the program produces cornerplots for 
%   the posterior distribution of parameter estimates. Default = 0.
%
%
% OUTPUT
%
% *BCIfitResults
% A structure with various fields containing fitting results from the BADS,
% VBMC and MCMC regimes:
% 
% - "data"
%   Structure with fields that contain the subject-specific data fields:
%   'trueLocsAV','responsesAV', and 'i_Conditions'.
%
% - "settings"
%   Structure containing all the settings (default and from user input)
%   used to run BCIfitModel.m. 
%
% - "BADS"
%   Structure containing the results of the Bayesian Adaptive Direct Search
%   optimization algorithm. It contains the following fields:
%   - "wrapper"   
%       A structure. By default BCIfitModel makes use of a wrapper function
%       around BADS that performs a Gridsearch and runs multiple BADS 
%       optimizations, to select the best result. This structure contains
%       information about the number of runs and the non-best solutions. 
%   - "fittedParams"
%       A vector of size [1 x nParams2Fit] with the parameter values as
%       returned by the best BADS solution (MAP or MLE depending on switch: 
%       "ForceLogLikelihood"). The order of the columns corresponds 
%       to the order of input parameters in "ParamNames2Fit".
%   - "exitflag"
%       A scalar with the exitflag of the best solution as returned by BADS
%       (its meaning is similar to exitflags of MATLAB's 'fminsearch').  
%   - "output"
%       The output structure of the best solution as returned by BADS 
%       containing more detailed information about the fitting procedure. 
%   - "prob"
%       A struct with three fields (scalars): 'logPrior', 'logLikelihood', 
%       and 'logPosterior'. These represent the probabilities associated 
%       with the BADS optimized parameter values. 
%   - "goodnessOfFit"
%       A structure with the following fields: 
%       (1) 'discreteCounts' a structure with arrays of response counts for
%       discrete responses (A, V, CSJ), and informational vectors of 
%       conditional info for the rows in each array (columns are the 
%       various response options: RespLocs or CSJ=1 / CSJ=2). These were 
%       used to compute the absolute goodness-of-fit absGoF (below). 
%       (2) 'predictedResps' a cell array (per condition) of structures
%       that contain the participant's true responses and BCI predicted 
%       response frequencies using the BADS optimized parameters. 
%       (3) 'LLvectors' a structure containing vectors of log-likehoods 
%       (per response) for the BCI model, the random-responses model, and 
%       the maximum-LL model (which simulates responses to maximise the LL
%       given the stochastic model structure). 
%       (4) 'InfoVectors' a structure containing vectors with trial 
%       information about the responses in LLvectors (i.e. stimulus 
%       locations, condition number, and participants' responses). 
%       (5) 'R2' a structure containing various R-squared measures of 
%       goodness-of-fit using the likelihood ratio-based method by 
%       Nagelkerke, 1991, Biometrika. An R2 of 0 (or negative) indicates 
%       that the BCI model performs no better (or worse!) than a model that
%       generates random responses, whereas an R2 of 1 indicates a perfect
%       fit (or as good as is possible given the stochastic nature of the
%       model). For an explanation of the various measures see the bottom
%       of "BCIcomputeGoF.m". (Use R2_Nagelkerke for continuous responses!)
%       (6) 'absGoF' a structure containing various absolute 
%       goodness-of-fit measures for discrete responses using the Kullback
%       Leibler divergence-based method by Acerbi, Dokka, Angelaki & Ma, 
%       2018, PLOS comp. Biology; see https://github.com/lacerbi/gofit. 
%       Again (same as R2 above), 0 means no better than chance, around 1
%       means a near perfect fit. (Use absGoF for discrete responses!)
%
% - "VBMC" [Optional, if VBMC was requested]
%   A structure with the following fields:
%   - "wrapper"   
%       A structure. By default BCIfitModel makes use of a wrapper function
%       around VBMC that runs multiple VBMC optimizations and selects the 
%       best result. This structure contains information about the number 
%       of runs and the non-best solutions. 
%   - "vp"
%       A structure containing the variational posterior solution (A 
%       mixture of Gaussians as an estimate of the posterior probability
%       distribution), but beware that this solution is in converted space
%       (e.g. sigmas are fitted in log-domain and values between 0 and 1 
%       are fitted in logit-domain; see 'BCIconvertParams2Fit.m')
%   - "elbo"
%       A scalar that represents the variational expected lower bound on
%       the log marginal likelihood (log model evidence). This measure can
%       be used for model comparison. See Gelman & Shalizi (British J. of 
%       Mathematical and Statistical Psychology, 2013) for the difference
%       between the ELBO (i.e. p(y)) and other (predictive) measures such
%       as PSISLOO (see MCMC below).
%   - "elbo_sd"
%       A scalar that represents the standard deviation of the elbo 
%       estimate (see above), as computed via Bayesian quadrature. This SD
%       should be small after convergence of VBMC.
%   - "exitflag"
%       A scalar with the exitflag as returned by the VBMC algorithm.
%   - "output"
%       A structure as returned by the VBMC algorithm.
%   - "validation"
%       A structure containing the results of a battery of diagnostics
%       tests that was run on all of the VBMC solutions. E.g. it compares
%       variational posterior solutions in terms of Kullback-Leibler 
%       divergence and total variation distance. If the exitflag of this
%       validation routine is 1 then multiple solutions converged to a
%       similar posterior distribution, thus making it more likely that
%       this solution is close to the real posterior distribution. 
%   - "stats" 
%       A structure containing the following descriptives of the posterior
%       distribution: 'post_mean' and 'post_std' (for each parameter), 
%       'post_iqr', a matrix with the 'quantiles' (rows) of each parameter
%       (column). Finally, 'post_cov' represents the covariance matrix of 
%       the posterior distribution over the parameters. But beware that 
%       these statistics are all in converted space (e.g. BCIfitModel 
%       sigmas are fitted in log-domain and values between 0 and 1 are 
%       fitted in logit-domain; see 'BCIconvertParams2Fit.m').
%
% - "MCMC" [Optional, if MCMC was requested]
%   A structure with the following fields:
%   - "mcmcParams"
%       A matrix of size [MCMCcount x nParams2Fit] with the parameter 
%       values of the MCMC procedure (order of params is defined by 
%       ParamNames2Fit).
%   - "prob"
%       A structure with 4 fields containing the computed probabilities 
%       for each of the MCMCcount samples: 'mcmcLogPrior', 
%       'mcmcLogLikelihood', 'logPosterior', and 'mcmcLL_responses'. 
%       The first three fields are vectors of size [MCMCcount x 1], whereas
%       the fourth field is a matrix of size [MCMCcount x nResponses], 
%       where nResponses can be larger than nTrials if more than one 
%       response is given per trial (e.g. localization and common source 
%       judgment). This matrix is used for computation of PSISLOO (see 
%       below). The probabilities in each row correspond to the parameter
%       values in mcmcParams.
%   - "converge"
%       A structure with three fields: (1) 'exitflag' tells you whether
%       MCMC has converged. (2) 'maxR' - The maximum (across parameters) 
%       of the Potential Scale Reduction Factor (see Gelman et al., 2013, 
%       Bayesian Data Analysis; chapter 11). (3) 'minNeff' - The minimum
%       (across parameters) of the Effective Number of Simulation Draws
%       (see Gelman et al., 2013, Bayesian Data Analysis; chapter 11).
%   - "output"
%       A structure containing information about the MCMC procedure: E.g.    
%       the overall number of function evaluations, the effective thinning
%       that was used, and convergence measures 'R' and 'Neff' split per
%       parameter.
%   - "goodnessOfFit"
%       A structure with the following fields: 
%       (1) 'PSISLOO', a structure containing the results from the Pareto
%       smoothed importance sampling leave-one-out function by Aki Vehtari 
%       https://github.com/avehtari/PSIS/blob/master/m/psisloo.m. 
%       E.g. PSISLOO.loo can be used for model comparison purposes (see: 
%       Aki Vehtari, Andrew Gelman and Jonah Gabry, 2017, Statistics and 
%       Computing). 
%       (2) 'R2' a structure containing various R-squared measures of 
%       goodness-of-fit using the method by Nagelkerke, 1991, Biometrika.
%       (Use R2_Nagelkerke for continuous responses)
%       (3) 'absGoF' a structure containing various absolute 
%       goodness-of-fit measures for discrete responses using the
%       method by Acerbi, Dokka, Angelaki & Ma, 2018, PLOS Comp. Biology; 
%       see https://github.com/lacerbi/gofit. 
%       (Use absGoF for discrete responses)
%       - N.B. for both (2) and (3) above we used the cross-validated 
%         log-likelihoods (i.e. PSISLOO estimates).     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check input arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add relevant subfolder to the path
me = mfilename;                                                 %what is my filename
pathstr = fileparts(which(me));                                 %get my location
addpath([pathstr filesep 'BCIscripts']);                        %add 'BCIscripts' folder to the path

%Is this a continuation from a previous run?
if nargin < 2
    assert(isstruct(trueLocsAV),'If only one input argument is provided, then it should be a BCIfitResults output structure from a previous run of BCIfitModel.m');
    BCIfitResults = trueLocsAV;
    
    if ~BCIfitResults.settings.BADS && ~isfield(BCIfitResults,'BADS')
        error('You first need to fit a maximum using BADS before calling the function separately for VBMC, MCMC or Display');
    end
    
%In case of a new run...     
else    
    [nTrials,nModalities] = size(trueLocsAV);
    assert(nModalities == 2,'trueLocsAV must have two columns: one for each sensory modality');

    [nResponses,nResponseTypes] = size(responsesAVC);
    assert(nTrials == nResponses,'responsesAVC must have the same number of rows as trueLocsAV: one for each trial');
    assert(nResponseTypes > 1,'responsesAVC must have at least two columns: one for each sensory modality');
    assert(nResponseTypes < 4,'responsesAVC cannot have more than three columns: one for each sensory modality and one for common-source judgments');
    if nResponseTypes == 2
        responsesAVC = [responsesAVC nan(nTrials,1)];
    end

    if nargin < 3
        i_Conditions = true(nTrials,1);
    elseif isempty(i_Conditions)
        i_Conditions = true(nTrials,1);
    else
        assert(nTrials == size(i_Conditions,1),'i_Conditions must have the same number of rows as trueLocsAV: one for each trial');
    end

    if nargin < 4
        OptionsStruct = struct([]);
    elseif isempty(OptionsStruct)
        OptionsStruct = struct([]);
    else
        assert(isstruct(OptionsStruct),'OptionsStruct must be a structure whose fields describe non-default settings');
    end
    
    %Evaluate 'OptionsStruct' and set default settings
    P = BCIsetOptions(OptionsStruct,responsesAVC);
    
    %Find all unique A,V location combinations for this subject
    P.sAsVconds = BCIfindLocCombinations(trueLocsAV);
    
    %Initialize output structure
    BCIfitResults = [];
    BCIfitResults.data.trueLocsAV = trueLocsAV;
    BCIfitResults.data.responsesAVC = responsesAVC;
    BCIfitResults.data.i_Conditions = i_Conditions;
    BCIfitResults.settings = P;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter optimization using BADS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.BADS
    BCIfitResults.BADS = BCIbads(BCIfitResults);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variational Bayesian Monte Carlo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.VBMC
    BCIfitResults.VBMC = BCIvbmc(BCIfitResults);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Markov Chain Monte Carlo %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.MCMC
    BCIfitResults.MCMC = BCImcmc(BCIfitResults);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Display Results %%%
%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.Display
    BCIdisplay(BCIfitResults);
end

end %[EOF]
