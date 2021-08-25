function BCIdisplayGroup(dataPath,savePath)
%Create summary figures for all 'BCIfitResults' in input argument dataPath.
%The folder identified by 'dataPath' should contain .mat files, where each
%.mat file contains a variable called "BCIfitResults" that is the output 
%structure of a previous call to BCIfitModel incl. BADS, VBMC, and/or MCMC.
%For each of these files the BCIdisplay.m function is called and the
%figures will be saved in 'savePath' (optional as 2nd input argument,
%otherwise 'savePath' = 'dataPath'). The savePath folder will be created if
%it does not exist yet. 

if nargin < 2
    savePath = dataPath;
else
    if 7~=exist(savePath,'dir')
       mkdir(savePath);
    end
end

%Ensure that BCIfitModel has been added to the path
me = mfilename;                                                             %what is my filename
myPath = fileparts(which(me));                                              %get my location
addpath(myPath);                                                            %add 'BCIscripts' to the path

%Gather the .mat filenames in dataPath
FileNames = dir(fullfile(dataPath,'*.mat'));
FileNames = {FileNames.name}';
nFiles = size(FileNames,1);

%Call display function for all files in a for-loop 
for iFile = 1:nFiles
    fileNameTmp = [dataPath filesep FileNames{iFile}];
    listOfVarsTmp = who('-file', fileNameTmp);
    if ismember('BCIfitResults', listOfVarsTmp) 
        load(fileNameTmp,'BCIfitResults');
        modifyBCIfitResults;                                                %Modify some settings in BCIfitResults (see nested helper function below) 
        displayAndSaveFigs(BCIfitResults,savePath,FileNames{iFile});        %See helper function below
        clear BCIfitResults                                                 %Free memory
    else
        disp(['File nr ' num2str(iFile) ' (' FileNames{iFile} ') does not contain a BCIfitResults variable. This file will be skipped.']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This nested function can access the workspace of the parent function
%It is memory efficient because there is no need to copy variables into input/output arguments.   
function modifyBCIfitResults

%Update settings to perform display only
BCIfitResults.settings.BADS = 0; 
BCIfitResults.settings.VBMC = 0; 
BCIfitResults.settings.MCMC = 0; 
BCIfitResults.settings.Display = 1; 

end %[End of Nested Function]

end %[End of Main/Parent Function]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%

%This function cannot access the workspace of the main function, but since it does not modify its input argument no additional memory is necessary   
function displayAndSaveFigs(BCIfitResults,savePath,filename)

figHandles = BCIdisplay(BCIfitResults);
nFigs = numel(figHandles);

%Get rid of '.mat' in filename
filename = filename(1:(end-4));

%Save and close all figures
for i=1:nFigs
    tmpName = get(figHandles(i),'Name');
    if tmpName(4) == '_'
        tmpName = tmpName(1:3);
    else
        tmpName = tmpName(1:4);
    end
    FullName = [savePath filesep filename '_FIG' num2str(i) '_' tmpName '.png'];
    saveas(figHandles(i),FullName);
    close(figHandles(i));
end

end %[EoF]

