%This little example script prepares a fake set of responses and then fits
%the BCI model to these responses. Feel free to play with the settings and
%recover the parameters yourself. Additionally, fitting the same dataset 
%multiple times may indicate how well the original parameters can be
%recovered. 

clear all;
close all;
clc

%Ensure that BCIfitModel has been added to the path
me = mfilename;                                                             %what is my filename
myPath = fileparts(which(me));                                              %get my location
BCIrootPathstr = fileparts(myPath);                                         %assume that this file is in a subfolder, e.g. BCIscripts            
addpath(BCIrootPathstr);                                                    %add 'BCIfitModel' to the path 
addpath([BCIrootPathstr filesep 'BCIscripts']);                             %add 'BCIscripts' subfolder to the path because we call "BCIgenerateResponsesMonteCarlo" directly (see below)

%Set the stimulus locations ('BCIgenerateDataset' will create a factorial design for sA,sV combinations that includes unisensory conditions)   
StimLocs = [-10 -5 0 5 10];
nTrialsPerAVcond = 100;

%Set the decision function, response locations and stimulus locations
P = cell(1);                                                                %P is a cell array of structs (one cell per experimental condition) where each struct indicates the non-default parameter settings
P{1}.decisionFun = 'ModelAveraging';                                        %'ModelAveraging', 'ModelSelection','ProbabilityMatching'
P{1}.RespLocs = StimLocs; %NaN;                                             %E.g. NaN for continuous responses or "StimLocs" for discrete responses

%Set some parameters to recover
P{1}.lapseR = 0.04;
P{1}.Pcommon = 0.33;
P{1}.sigmaP = 24;
P{1}.sigmaA = 7.5;
P{1}.sigmaV = 2;
%P{1}.deltaSigmaV = 0;
%P{1}.deltaXA = 0.2;
%P{1}.deltaXV = -0.1;

%Add a second condition wherein the visual reliability has decreased
P{2} = P{1};                                                               
P{2}.sigmaV = 5;

%Generate a dataset for one participant
[trueLocsAV,responsesAVC,i_Conditions] = BCIgenerateDataset(P,StimLocs,nTrialsPerAVcond);
responsesAVC(:,3) = [];                                                     %Optionally delete the common source judgments

%Prepare some settings for BCIfitModel
OptionsStruct = [];
OptionsStruct.Display = 1;                                                  %Plot some results at the end
OptionsStruct.VBMC = 1; OptionsStruct.MCMC = 1;                             %Besides BADS optimization (default) also run VBMC and MCMC        
OptionsStruct.nAttempts.BADS = [4 8]; OptionsStruct.nAttempts.VBMC = [4 8]; OptionsStruct.nAttempts.MCMC = [1 1];   %[Minimum number of attempts, maximum number of attempts] - after which the best result is selected
OptionsStruct.parallel.BADS = 1; OptionsStruct.parallel.VBMC = 1; OptionsStruct.parallel.MCMC = 0;                  %Use Matlab Parallelization toolbox? (not recommended for MCMC)

OptionsStruct.decisionFun = P{1}.decisionFun;
OptionsStruct.RespLocs = P{1}.RespLocs; 

OptionsStruct.integrateMethod = 2;                                          %1 or 2, for continuous responses only. Use mean motor noise pdf or an estimate of the actual pdf (the second is faster but slightly less precise). 
OptionsStruct.triangpdfSpeedUp = 0;                                         %Use a triangular approximation of normal pdf for motor noise for a small speed up (only relevant for integration method 1)

%OptionsStruct.ParamNames2Fit = {'lapseR','Pcommon','sigmaP','sigmaA','sigmaV','deltaSigmaV','deltaXA','deltaXV'}; %<-- feel free to add others (or add a second condition as below)

OptionsStruct.ParamNames2Fit = {'lapseR','Pcommon','sigmaP','sigmaA','sigmaV','sigmaV'};   %Note that sigmaV appears twice (one for each experimental condition)
OptionsStruct.ParamsPerCond = {[1 2 3 4 5],[1 2 3 4 6]};                                   %The first sigmaV belongs to the first condition, the second sigmaV belongs to the second condition. 
                                                                                           %Other parameters are shared across conditions
%Finally, call the BCI model fit function (this may take a few minutes!)
BCIfitResults = BCIfitModel(trueLocsAV,responsesAVC,i_Conditions,OptionsStruct);            
                                                                                            
%Please do take a look at the returned structure 'BCIfitResults' !!!
