function [BCIfitResultsGroup,FileNames,BCIfitResultsGroupMean,fig_handles] = BCIgatherGroupSummary(dataPath,GroupMeanBADSplotBool)
%Create one datafile with summary statistics for all 'BCIfitResults' in the
%folder that is identified by input argument 'dataPath'. This folder should
%contain .mat files, where each .mat file contains a variable called 
%"BCIfitResults" that is the output structure of a previous call to 
%BCIfitModel incl. BADS, VBMC, and/or MCMC. This function assumes that each
%of these files was created using the same settings (but different data). 
%For each of these files the fitted parameters etc. will be collected. 
%The collected results are returned as one summary struct. The second 
%output contains all FileNames and they are listed in the same order as the
%summary results (alphabetical according to Matlab's 'dir'). The third
%output is a BCIfitResults struct with the group mean results. The fourth
%output contains the figure handles to the figures produced by BCIdisplay.

if nargin < 2
    GroupMeanBADSplotBool = 1;
end
if ~GroupMeanBADSplotBool
    fig_handles = [];       %Ensure that output parameter exists
end

%Ensure that BCIfitModel has been added to the path
me = mfilename;                                                             %what is my filename
myPath = fileparts(which(me));                                              %get my location
addpath(myPath);                                                            %add 'BCIscripts' to the path

%Gather the .mat filenames in dataPath
FileNames = dir(fullfile(dataPath,'*.mat'));
FileNames = {FileNames.name}';
nFiles = size(FileNames,1);

%Find all files that contain a BCIfitResults structure
useFile = true(nFiles,1);
for iFile=1:nFiles
    fileNameTmp = [dataPath filesep FileNames{iFile}];
    listOfVarsTmp = who('-file', fileNameTmp);
    if ~ismember('BCIfitResults', listOfVarsTmp) 
        useFile(iFile) = false;
        disp(['File nr ' num2str(iFile) ' (' FileNames{iFile} ') does not contain a BCIfitResults variable. This file will be skipped.']);
    end
end
nFiles = sum(useFile);
FileNames = FileNames(useFile,1);
if nFiles == 0; error('There were no files that contained a BCIfitResults structure in dataPath'); end

%Initialize the summary structures
BCIfitResultsGroup.nSubjects = nFiles;
BCIfitResultsGroup.settings = [];
BCIfitResultsGroup.BADS = [];
BCIfitResultsGroup.VBMC = [];
BCIfitResultsGroup.MCMC = [];

%Call the collect functions for all files in a for-loop 
for iFile = 1:nFiles
    fileNameTmp = [dataPath filesep FileNames{iFile}];
    load(fileNameTmp,'BCIfitResults');
    BCIfitResultsGroup.settings = collectSettings(BCIfitResultsGroup.settings,BCIfitResults.settings,iFile,FileNames);
    if isfield(BCIfitResults,'BADS')
        BCIfitResultsGroup.BADS = collectBADS(BCIfitResultsGroup.BADS,BCIfitResults.BADS,iFile,nFiles);
    end
    if isfield(BCIfitResults,'VBMC')
        BCIfitResultsGroup.VBMC = collectVBMC(BCIfitResultsGroup.VBMC,BCIfitResults.VBMC,iFile,nFiles);
    end
    if isfield(BCIfitResults,'MCMC')
        BCIfitResultsGroup.MCMC = collectMCMC(BCIfitResultsGroup.MCMC,BCIfitResults.MCMC,iFile,nFiles);
    end
    disp(['File #' num2str(iFile) ' out of ' num2str(nFiles) ' done..']);
end

%Delete the structure that were not filled
if isempty(BCIfitResultsGroup.BADS)
    BCIfitResultsGroup = rmfield(BCIfitResultsGroup,'BADS');
end
if isempty(BCIfitResultsGroup.VBMC)
    BCIfitResultsGroup = rmfield(BCIfitResultsGroup,'VBMC');
end
if isempty(BCIfitResultsGroup.MCMC)
    BCIfitResultsGroup = rmfield(BCIfitResultsGroup,'MCMC');
end

%Create one group-mean BCIfitResults structure with BADS only and plot these results   
if isfield(BCIfitResults,'BADS')
    BCIfitResultsGroupMean = compGroupMeanBADS(BCIfitResultsGroup);
end

%Plot the group mean (for BADS only)
if GroupMeanBADSplotBool
    fig_handles = BCIdisplay(BCIfitResultsGroupMean);
end

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

function groupSettings = collectSettings(groupSettings,subjSettings,iSubj,FileNames)
    
    %Move fixed parameters into separate group 'fixedParams' and delete NaN entries (this step is unnecessary because we use isequaln instead of isequal below...)   
    if ismember('lapseR',subjSettings.ParamNames2Fit)                                                                           
        fixedParamNamesAll = {'Pcommon','sigmaA','sigmaV','sigmaP','muP','deltaSigmaA','deltaSigmaV','deltaXA','deltaXV','CSJthresh','sigmaMotor'};             %lapseR is fit and overwrites lapseRA, lapseRV and lapseRC
        subjSettings = rmfield(subjSettings,{'lapseR','lapseRA','lapseRV','lapseRC'});
    elseif ~isnan(subjSettings.lapseR)                                                                                          
        fixedParamNamesAll = {'Pcommon','sigmaA','sigmaV','sigmaP','muP','deltaSigmaA','deltaSigmaV','deltaXA','deltaXV','CSJthresh','sigmaMotor','lapseR'};    %lapseR is fixed and overwrites lapseRA, lapseRV and lapseRC
        subjSettings = rmfield(subjSettings,{'lapseRA','lapseRV','lapseRC'});
    else                                                                                                                     %lapseR is undefined (NaN) and does not overwrite fixed values for lapseRA, lapseRV and lapseRC
        fixedParamNamesAll = {'Pcommon','sigmaA','sigmaV','sigmaP','muP','deltaSigmaA','deltaSigmaV','deltaXA','deltaXV','CSJthresh','sigmaMotor','lapseRA','lapseRV','lapseRC'};
        subjSettings = rmfield(subjSettings,'lapseR');
    end
    for i=1:length(fixedParamNamesAll) 
        if ~ismember(fixedParamNamesAll{i},subjSettings.ParamNames2Fit)
            subjSettings.fixedParams.(fixedParamNamesAll{i}) = subjSettings.(fixedParamNamesAll{i});
        end
        subjSettings = rmfield(subjSettings,fixedParamNamesAll{i});
    end
    fnamesTmp = fieldnames(subjSettings);                                   
    subjSettings = orderfields(subjSettings,[1:4 length(fnamesTmp) 5:(length(fnamesTmp)-1)]);   %Reorder the fieldnames such that 'fixedParams' (added as last) gets place number 5
    
    %First subject: simple copy
    if isempty(groupSettings)
        groupSettings = subjSettings;
        
    %For all other subjects: assure that the settings are equal    
    else
        if ~isequaln(groupSettings,subjSettings) %n.b. isequaln considers NaNs equal to each other
            error(['Settings for file #' num2str(iSubj) ' (' FileNames{iSubj} ') are not equal to the settings of file #1 (' FileNames{1} '). Please ensure that all files were run with identical options for BCIfitModel']);
        end
    end
end

%-------------------------------------------------------------------------

function groupBADS = collectBADS(groupBADS,subjBADS,iSubj,nSubjects)

    %First subject: initialize everything
    if isempty(groupBADS)
        groupBADS.wrapper.bestIdx = nan(nSubjects,1);
        groupBADS.wrapper.fittedParams = nan([size(subjBADS.wrapper.fittedParams) nSubjects]);
        groupBADS.wrapper.fittedValue = nan([length(subjBADS.wrapper.fittedValue) nSubjects]);
        groupBADS.wrapper.exitflags = nan([length(subjBADS.wrapper.fittedValue) nSubjects]);
        
        groupBADS.fittedParams = nan(nSubjects,numel(subjBADS.fittedParams));
        groupBADS.exitflag = nan(nSubjects,1);
        
        groupBADS.prob.logPrior = nan(nSubjects,1);
        groupBADS.prob.logLikelihood = nan(nSubjects,1);
        groupBADS.prob.logPosterior = nan(nSubjects,1);
        
        groupBADS.goodnessOfFit.discreteCounts = cell(nSubjects,1);
        groupBADS.goodnessOfFit.predictedResps = cell([size(subjBADS.goodnessOfFit.predictedResps) nSubjects]);
        groupBADS.goodnessOfFit.R2 = cell(nSubjects+1,size(subjBADS.goodnessOfFit.R2,2)); 
        groupBADS.goodnessOfFit.R2(1,:) = subjBADS.goodnessOfFit.R2(1,:);
        groupBADS.goodnessOfFit.absGoF = cell(nSubjects+1,size(subjBADS.goodnessOfFit.absGoF,2));
        groupBADS.goodnessOfFit.absGoF(1,:) = subjBADS.goodnessOfFit.absGoF(1,:);
        
        groupBADS.goodnessOfFit.nResponses = nan(nSubjects,1);
    end
        
    %All subjects: Fill the data with this subject
    groupBADS.wrapper.bestIdx(iSubj,1) = subjBADS.wrapper.bestIdx;
    groupBADS.wrapper.fittedParams(:,:,iSubj) = subjBADS.wrapper.fittedParams;
    groupBADS.wrapper.fittedValue(:,iSubj) = subjBADS.wrapper.fittedValue;
    if isfield(subjBADS.wrapper,'allResults')
        groupBADS.wrapper.exitflags(:,iSubj) = subjBADS.wrapper.allResults.exitflags;            %Compatibility issues
    elseif isfield(subjBADS.wrapper,'exitflags')
        groupBADS.wrapper.exitflags(:,iSubj) = subjBADS.wrapper.exitflags;
    end
    
    groupBADS.fittedParams(iSubj,:) = subjBADS.fittedParams;
    groupBADS.exitflag(iSubj,1) = subjBADS.exitflag;

    groupBADS.prob.logPrior(iSubj,1) = subjBADS.prob.logPrior;
    groupBADS.prob.logLikelihood(iSubj,1) = subjBADS.prob.logLikelihood;
    groupBADS.prob.logPosterior(iSubj,1) = subjBADS.prob.logPosterior;
    
    groupBADS.goodnessOfFit.discreteCounts{iSubj,1} = subjBADS.goodnessOfFit.discreteCounts;
    groupBADS.goodnessOfFit.predictedResps(:,:,iSubj) = subjBADS.goodnessOfFit.predictedResps;
    groupBADS.goodnessOfFit.R2(iSubj+1,:) = subjBADS.goodnessOfFit.R2(2,:);
    groupBADS.goodnessOfFit.absGoF(iSubj+1,:) = subjBADS.goodnessOfFit.absGoF(2,:);
    
    %Gather total number of responses in order to compute BIC (Bayesian information criterion) 
    groupBADS.goodnessOfFit.nResponses(iSubj,1) = sum(~isnan(subjBADS.goodnessOfFit.InfoVectors.Aresps))+sum(~isnan(subjBADS.goodnessOfFit.InfoVectors.Vresps))+sum(~isnan(subjBADS.goodnessOfFit.InfoVectors.CSJresps));
    
end %[EoF]

%-------------------------------------------------------------------------

function groupVBMC = collectVBMC(groupVBMC,subjVBMC,iSubj,nSubjects)

    %First subject: initialize everything
    if isempty(groupVBMC)
        groupVBMC.wrapper.idx_best = nan(nSubjects,1);
        groupVBMC.wrapper.paramMeans = nan([size(subjVBMC.wrapper.paramMeans) nSubjects]);
        groupVBMC.wrapper.paramSDs = nan([size(subjVBMC.wrapper.paramSDs) nSubjects]);
        groupVBMC.wrapper.elbo = nan([length(subjVBMC.wrapper.elbo) nSubjects]);
        groupVBMC.wrapper.elbo_sd = nan([length(subjVBMC.wrapper.elbo) nSubjects]);
        groupVBMC.wrapper.exitflag = nan([length(subjVBMC.wrapper.elbo) nSubjects]);
        
        groupVBMC.elbo = nan(nSubjects,1);
        groupVBMC.elbo_sd = nan(nSubjects,1);
        groupVBMC.exitflag = nan(nSubjects,1);
        
        groupVBMC.validation.sKL_best = nan([length(subjVBMC.validation.sKL_best) nSubjects]);
        groupVBMC.validation.maxmtv_best = nan([length(subjVBMC.validation.maxmtv_best) nSubjects]);
        groupVBMC.validation.exitflag = nan(nSubjects,1);
        
        groupVBMC.stats.quantiles = subjVBMC.stats.quantiles;
        groupVBMC.stats.post_iqr = nan([size(subjVBMC.stats.post_iqr) nSubjects]);
        groupVBMC.stats.post_mean = nan([length(subjVBMC.stats.post_mean) nSubjects]);
        groupVBMC.stats.post_std = nan([length(subjVBMC.stats.post_std) nSubjects]);
        groupVBMC.stats.post_cov = nan([size(subjVBMC.stats.post_cov) nSubjects]);
    end 
        
    %All subjects: Fill the data with this subject
    groupVBMC.wrapper.idx_best(iSubj,1) = subjVBMC.wrapper.idx_best;
    groupVBMC.wrapper.paramMeans(:,:,iSubj) = subjVBMC.wrapper.paramMeans;
    groupVBMC.wrapper.paramSDs(:,:,iSubj) = subjVBMC.wrapper.paramSDs;
    groupVBMC.wrapper.elbo(:,iSubj) = subjVBMC.wrapper.elbo;
    groupVBMC.wrapper.elbo_sd(:,iSubj) = subjVBMC.wrapper.elbo_sd;
    groupVBMC.wrapper.exitflag(:,iSubj) = subjVBMC.wrapper.exitflag;
    
    groupVBMC.elbo(iSubj,1) = subjVBMC.elbo;
    groupVBMC.elbo_sd(iSubj,1) = subjVBMC.elbo_sd;
    groupVBMC.exitflag(iSubj,1) = subjVBMC.exitflag;

    groupVBMC.validation.sKL_best(:,iSubj) = subjVBMC.validation.sKL_best;
    groupVBMC.validation.maxmtv_best(:,iSubj) = subjVBMC.validation.maxmtv_best;
    groupVBMC.validation.exitflag(iSubj,1) = subjVBMC.validation.exitflag;
        
    groupVBMC.stats.post_iqr(:,:,iSubj) = subjVBMC.stats.post_iqr;
    groupVBMC.stats.post_mean(:,iSubj) = subjVBMC.stats.post_mean;
    groupVBMC.stats.post_std(:,iSubj) = subjVBMC.stats.post_std;
    groupVBMC.stats.post_cov(:,:,iSubj) = subjVBMC.stats.post_cov;
    
end %[EoF]

%-------------------------------------------------------------------------

function groupMCMC = collectMCMC(groupMCMC,subjMCMC,iSubj,nSubjects)

    %First subject: initialize everything
    if isempty(groupMCMC)
        groupMCMC.mcmcParams = nan([size(subjMCMC.mcmcParams) nSubjects]);
        
        groupMCMC.prob.mcmcLogPrior = nan([length(subjMCMC.prob.mcmcLogPrior) nSubjects]);
        groupMCMC.prob.mcmcLogLikelihood = nan([length(subjMCMC.prob.mcmcLogLikelihood) nSubjects]);
        
        groupMCMC.converge.exitflag = nan(nSubjects,1);
        groupMCMC.converge.maxR = nan(nSubjects,1);
        groupMCMC.converge.minNeff = nan(nSubjects,1);
        
        groupMCMC.output.thin = nan(nSubjects,1);
        groupMCMC.output.funccount = nan(nSubjects,1);
        groupMCMC.output.nslicecollapsed = nan(nSubjects,1);
        groupMCMC.output.R = nan([length(subjMCMC.output.R) nSubjects]);
        groupMCMC.output.Neff = nan([length(subjMCMC.output.Neff) nSubjects]);
        groupMCMC.output.tau = nan([length(subjMCMC.output.tau) nSubjects]);
        
        groupMCMC.goodnessOfFit.PSISLOO.loo = nan(nSubjects,1);
        groupMCMC.goodnessOfFit.PSISLOO.kn0 = nan(nSubjects,1);
        groupMCMC.goodnessOfFit.PSISLOO.kn1 = nan(nSubjects,1);
        groupMCMC.goodnessOfFit.PSISLOO.kn2 = nan(nSubjects,1);
        
        groupMCMC.goodnessOfFit.R2 = cell(nSubjects+1,size(subjMCMC.goodnessOfFit.R2,2)); 
        groupMCMC.goodnessOfFit.R2(1,:) = subjMCMC.goodnessOfFit.R2(1,:);
        groupMCMC.goodnessOfFit.absGoF = cell(nSubjects+1,size(subjMCMC.goodnessOfFit.absGoF,2));
        groupMCMC.goodnessOfFit.absGoF(1,:) = subjMCMC.goodnessOfFit.absGoF(1,:);
    end 
        
    %All subjects: Fill the data with this subject
    groupMCMC.mcmcParams(:,:,iSubj) = subjMCMC.mcmcParams;
        
    groupMCMC.prob.mcmcLogPrior(:,iSubj) = subjMCMC.prob.mcmcLogPrior;
    groupMCMC.prob.mcmcLogLikelihood(:,iSubj) = subjMCMC.prob.mcmcLogLikelihood;

    groupMCMC.converge.exitflag(iSubj,1) = subjMCMC.converge.exitflag;
    groupMCMC.converge.maxR(iSubj,1) = subjMCMC.converge.maxR;
    groupMCMC.converge.minNeff(iSubj,1) = subjMCMC.converge.minNeff;

    groupMCMC.output.thin(iSubj,1) = subjMCMC.output.thin;
    groupMCMC.output.funccount(iSubj,1) = subjMCMC.output.funccount(2);
    groupMCMC.output.nslicecollapsed(iSubj,1) = subjMCMC.output.nslicecollapsed;
    groupMCMC.output.R(:,iSubj) = subjMCMC.output.R;
    groupMCMC.output.Neff(:,iSubj) = subjMCMC.output.Neff;
    groupMCMC.output.tau(:,iSubj) = subjMCMC.output.tau;

    groupMCMC.goodnessOfFit.PSISLOO.loo(iSubj,1) = subjMCMC.goodnessOfFit.PSISLOO.loo;
    groupMCMC.goodnessOfFit.PSISLOO.kn0(iSubj,1) = subjMCMC.goodnessOfFit.PSISLOO.kn0;
    groupMCMC.goodnessOfFit.PSISLOO.kn1(iSubj,1) = subjMCMC.goodnessOfFit.PSISLOO.kn1;
    groupMCMC.goodnessOfFit.PSISLOO.kn2(iSubj,1) = subjMCMC.goodnessOfFit.PSISLOO.kn2;
    
    groupMCMC.goodnessOfFit.R2(iSubj+1,:) = subjMCMC.goodnessOfFit.R2(2,:);
    groupMCMC.goodnessOfFit.absGoF(iSubj+1,:) = subjMCMC.goodnessOfFit.absGoF(2,:);
    
end %[EoF]

%-------------------------------------------------------------------------

function groupMeanResults = compGroupMeanBADS(BCIfitResultsGroup)
    
    %Copy settings
    groupMeanResults.nSubjects = BCIfitResultsGroup.nSubjects;
    groupMeanResults.settings = BCIfitResultsGroup.settings;
    discreteBool = all(~isnan(groupMeanResults.settings.RespLocs));
    respRange = groupMeanResults.settings.RespRange;
    
    %Compute some simple means for the BADS summary measures
    groupMeanResults.BADS.fittedParams = mean(BCIfitResultsGroup.BADS.fittedParams);
    groupMeanResults.BADS.exitflag = min(BCIfitResultsGroup.BADS.exitflag);                 %<-- minimum exitflag (not mean)

    groupMeanResults.BADS.prob.logPrior = mean(BCIfitResultsGroup.BADS.prob.logPrior);
    groupMeanResults.BADS.prob.logLikelihood = mean(BCIfitResultsGroup.BADS.prob.logLikelihood);
    groupMeanResults.BADS.prob.logPosterior = mean(BCIfitResultsGroup.BADS.prob.logPosterior);

    %Compute mean R2 and absGoF
    nConds = size(BCIfitResultsGroup.BADS.goodnessOfFit.R2,2);
    groupMeanResults.BADS.goodnessOfFit.R2 = cell(2,nConds);
    groupMeanResults.BADS.goodnessOfFit.R2(1,:) = BCIfitResultsGroup.BADS.goodnessOfFit.R2(1,:);
    groupMeanResults.BADS.goodnessOfFit.absGoF = cell(2,nConds);
    groupMeanResults.BADS.goodnessOfFit.absGoF(1,:) = BCIfitResultsGroup.BADS.goodnessOfFit.absGoF(1,:);
    for i=1:nConds
        groupMeanResults.BADS.goodnessOfFit.R2{2,i} = meanOfEachField(BCIfitResultsGroup.BADS.goodnessOfFit.R2(2:end,i));
        groupMeanResults.BADS.goodnessOfFit.absGoF{2,i} = meanOfEachField(BCIfitResultsGroup.BADS.goodnessOfFit.absGoF(2:end,i));
    end
    
    %Deal with predictedResps
    groupMeanResults.BADS.goodnessOfFit.predictedResps = meanPredictedResps(BCIfitResultsGroup.BADS.goodnessOfFit.predictedResps,discreteBool,respRange);
    
    %Don't compute means for VBMC and MCMC here.. This function is 'only' for being able to plot mean fitting results
    groupMeanResults.settings.BADS = 1; groupMeanResults.settings.VBMC = 0; groupMeanResults.settings.MCMC = 0; 
    
    %However, do include an analysis of the mean R2 and absGoF (since we have the little "meanOfEachField" function here)   
    if isfield(BCIfitResultsGroup,'MCMC')
        groupMeanResults.MCMC.goodnessOfFit.R2 = cell(2,nConds);
        groupMeanResults.MCMC.goodnessOfFit.R2(1,:) = BCIfitResultsGroup.MCMC.goodnessOfFit.R2(1,:);
        groupMeanResults.MCMC.goodnessOfFit.absGoF = cell(2,nConds);
        groupMeanResults.MCMC.goodnessOfFit.absGoF(1,:) = BCIfitResultsGroup.MCMC.goodnessOfFit.absGoF(1,:);
        for i=1:nConds
            groupMeanResults.MCMC.goodnessOfFit.R2{2,i} = meanOfEachField(BCIfitResultsGroup.MCMC.goodnessOfFit.R2(2:end,i));
            groupMeanResults.MCMC.goodnessOfFit.absGoF{2,i} = meanOfEachField(BCIfitResultsGroup.MCMC.goodnessOfFit.absGoF(2:end,i));
        end
    end
    
end %[EoF]

%-------------------------------------------------------------------------

function meanStruct = meanOfEachField(cellArrayWithStructs,meanStruct,fieldNameCell)

%First call
if nargin < 3
    meanStruct = cellArrayWithStructs{1};
    fieldNamesTmp = fieldnames(meanStruct);
    for i=1:length(fieldNamesTmp)
        meanStruct = meanOfEachField(cellArrayWithStructs,meanStruct,fieldNamesTmp(i));     
    end
else
    %Recursively call one-self for each structure field
    if isstruct(getfield(meanStruct,fieldNameCell{:}))
        fieldNamesTmp = fieldnames(getfield(meanStruct,fieldNameCell{:}));
        for i=1:length(fieldNamesTmp)
            meanStruct = meanOfEachField(cellArrayWithStructs,meanStruct,[fieldNameCell; fieldNamesTmp(i)]);     
        end
    %Evaluate the field: compute mean of all values across all structures in cell array      
    else
        nCells = numel(cellArrayWithStructs);
        collect = nan(1,nCells);
        for i=1:nCells
            collect(i) = getfield(cellArrayWithStructs{i},fieldNameCell{:});
        end
        meanStruct = setfield(meanStruct,fieldNameCell{:},mean(collect));   %simple mean across collected vector values
    end
end

end %[EoF]

%-------------------------------------------------------------------------

function meanCell = meanPredictedResps(cellArrayWithStructs,discreteBool,respRange)

%Initialize
[nAVlocConds,nConds,nSubjects] = size(cellArrayWithStructs);
meanCell = cellArrayWithStructs(:,:,1);

%Loop through all conditions and concatenate the responses of all subjects
for i=1:nAVlocConds
    for j=1:nConds
        for k=2:nSubjects
            meanCell{i,j}.respA = [meanCell{i,j}.respA; cellArrayWithStructs{i,j,k}.respA];
            meanCell{i,j}.respV = [meanCell{i,j}.respV; cellArrayWithStructs{i,j,k}.respV];
            meanCell{i,j}.respCSJ = [meanCell{i,j}.respCSJ; cellArrayWithStructs{i,j,k}.respCSJ];
            
            %Average the likelihood functions in case of discrete responses
            if discreteBool
                meanCell{i,j}.likeFunA = meanCell{i,j}.likeFunA + cellArrayWithStructs{i,j,k}.likeFunA;         %First sum..
                meanCell{i,j}.likeFunV = meanCell{i,j}.likeFunV + cellArrayWithStructs{i,j,k}.likeFunV;
            end
            meanCell{i,j}.likeFunCSJ = meanCell{i,j}.likeFunCSJ + cellArrayWithStructs{i,j,k}.likeFunCSJ;
        end
        if discreteBool
            meanCell{i,j}.likeFunA = meanCell{i,j}.likeFunA / nSubjects;                                        %..then normalize      
            meanCell{i,j}.likeFunV = meanCell{i,j}.likeFunV / nSubjects;
        end
        meanCell{i,j}.likeFunCSJ = meanCell{i,j}.likeFunCSJ / nSubjects;
    end
end

%Loop through all conditions and average the likelihood functions
if ~discreteBool
    nGrid = 1000;
    grid = linspace(respRange(1),respRange(2),nGrid);
    deltaX = (grid(nGrid)-grid(1))/(nGrid-1);
    for i=1:nAVlocConds
        for j=1:nConds
            meanCell{i,j}.likeFunA_grid = grid;
            meanCell{i,j}.likeFunV_grid = grid;
            
            meanCell{i,j}.likeFunA = lininterp1(cellArrayWithStructs{i,j,1}.likeFunA_grid,cellArrayWithStructs{i,j,k}.likeFunA,grid,[],deltaX);
            meanCell{i,j}.likeFunV = lininterp1(cellArrayWithStructs{i,j,1}.likeFunV_grid,cellArrayWithStructs{i,j,k}.likeFunV,grid,[],deltaX);
            
            for k=2:nSubjects
                meanCell{i,j}.likeFunA = meanCell{i,j}.likeFunA + lininterp1(cellArrayWithStructs{i,j,1}.likeFunA_grid,cellArrayWithStructs{i,j,k}.likeFunA,grid,[],deltaX);        %First sum ..
                meanCell{i,j}.likeFunV = meanCell{i,j}.likeFunV + lininterp1(cellArrayWithStructs{i,j,1}.likeFunV_grid,cellArrayWithStructs{i,j,k}.likeFunV,grid,[],deltaX);
            end
            meanCell{i,j}.likeFunA = meanCell{i,j}.likeFunA / nSubjects;                                                                                                            %.. then normalize      
            meanCell{i,j}.likeFunV = meanCell{i,j}.likeFunV / nSubjects;
        end
    end
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
if isnan(extrap)
    flag3 = (Xq < X(1)) | (Xq > X(end)) | isnan(Xq);
else
    flag3 = isnan(Xq);
    if isempty(extrap)
        flag1 = Xq < X(1);
        flag2 = Xq > X(end);
    else
        flag12 = (Xq < X(1)) | (Xq > X(end));
    end
end

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
out1 = V(Xqi);
out2 = V(Xqi+1);
Vout = (1-delta).*out1 + delta.*out2;

%Set values outside of X (extrapolation)
if ~isnan(extrap)
    %If extrap is set to empty: Xq outside X get V values of nearest edge.
    if isempty(extrap)                      
        Vout(flag1) = V(1);
        Vout(flag2) = V(nGrid);
    %If extrap is set to any scalar: Xq outside X get that specified value.
    else
        Vout(flag12) = extrap;
    end
end
Vout(flag3) = NaN; 

end %[EoF]
