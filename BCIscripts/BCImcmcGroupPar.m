function BCImcmcGroupPar(dataPath,maxAttemptsPerFile,appendBool)
%Run MCMC for multiple subjects in parallel. The folder identified by input
%argument dataPath should contain .mat files, where each .mat file contains
%a variable called "BCIfitResults" that is the output structure of a
%previous call to BCIfitModel incl. BADS and/or VBMC. These .mat files will
%be -appended with the structure field MCMC added to BCIfitResults.
%Also see the nested helper function 'modifyBCIfitResults' to make changes
%to the settings of the individual BCIfitResults files.

%Set a maximal number of attempts per file (files that don't converge at once will receive more attempts up to this maximum, or potentially even more as long as other files are still running too)
if nargin < 2
    maxAttemptsPerFile = 1;
end

%Should we -append the BCIfitResults output structure or save as '-v7.3' 
if nargin < 3
    appendBool = 0;                                                         %I.e. by default we simply save as v7.3 (ensuring compression but risk overwriting other variables)    
end

%Ensure that BCIfitModel has been added to the path
me = mfilename;                                                             %what is my filename
myPath = fileparts(which(me));                                              %get my location
addpath(fileparts(myPath));                                                 %add 'BCIfitModel' to the path (assuming that this file is in a subfolder, e.g. BCIscripts)

%Set up a parallel pool (in case it hasn't yet)
nParallelWorkers = feature('numcores');
currPool = gcp('nocreate');
if isempty(currPool)
    try %Handle occasional errors..                                         %"Error using parpool (line 145)
        currPool = parpool('local',nParallelWorkers);                       %Parallel pool failed to start with the following error. For more detailed
    catch %Try again..                                                      %information, validate the profile 'local' in the Cluster Profile Manager.
        pause(5);                                                           %Error using parallel.internal.pool.InteractiveClient>iThrowWithCause (line 670)
        myCluster = parcluster('local');                                    %Failed to initialize the interactive session. Error using  
        delete(myCluster.Jobs);                                             %parallel.internal.pool.InteractiveClient>iThrowIfBadParallelJobStatus (line 810)    
        currPool = parpool('local',nParallelWorkers);                       %The interactive communicating job finished with no message."  
    end                                                                     
end                                                                         %Also: Remove all jobs created with profile local (these are saved after a crash)

%Shuffle the random number generator once
rng('shuffle');

%Gather the .mat filenames in dataPath
FileNames = dir(fullfile(dataPath,'*.mat'));
FileNames = {FileNames.name}';
nFiles = size(FileNames,1);

%Initialize some trackers for each .mat file
fileConverged = zeros(1,nFiles);
fileCounts = zeros(2,nFiles);                                               %first row is for nCalled, second row is for nFinished
fileMaxRs = nan(1,nFiles);   

%Set an arbitrary large number of maximum total number of runs (only used for initializing the parallel tasks below)   
maxRuns = 2*maxAttemptsPerFile*nFiles+nParallelWorkers; 

%Initialize all tasks
FuturesArray(1,maxRuns) = parallel.FevalFuture();                         
fileNrsRunningNow = nan(1,maxRuns);
NotStartedBools = zeros(1,maxRuns);
StartTimes = zeros(1,maxRuns,'uint64');

%Randomly pick the first batch of fileNrs to run
shuffled_fileNrs_tooMany = repmat(randperm(nFiles),[1 ceil(nParallelWorkers/nFiles)]);      %replicate enough (or too many) of the shuffled file nrs
fileNrsRunningNow(1:nParallelWorkers) = shuffled_fileNrs_tooMany(:,1:nParallelWorkers);     %make sure there are not too many anymore

%Start MCMC on all parallel workers 
for newIdx = 1:nParallelWorkers
    StartTimes(newIdx) = tic;
    fileCounts(1,fileNrsRunningNow(newIdx)) = fileCounts(1,fileNrsRunningNow(newIdx))+1;
    fileNameTmp = [dataPath filesep FileNames{fileNrsRunningNow(newIdx)}];
    listOfVarsTmp = who('-file', fileNameTmp);
    if ismember('BCIfitResults', listOfVarsTmp) 
        load(fileNameTmp,'BCIfitResults');
        if convergedAlready(BCIfitResults)
            if ~fileConverged(fileNrsRunningNow(newIdx))                    %If this is the first time that this file was checked (it may have been chosen multiple times, but we want the message only once)
                disp(['File nr ' num2str(fileNrsRunningNow(newIdx)) ' (' FileNames{fileNrsRunningNow(newIdx)} ') has already converged. This file will be skipped.']);
                fileConverged(fileNrsRunningNow(newIdx)) = 1;
                fileCounts(:,fileNrsRunningNow(newIdx)) = inf;              %Ensure that this file is not chosen to run again
            end
            NotStartedBools(newIdx) = 1;
        else
            modifyBCIfitResults;                                            %Modify some settings in BCIfitResults (see nested helper function below) 
            FuturesArray(newIdx) = parfeval(@BCIfitModel,1,BCIfitResults);  %Expect 1 output argument, and use 1 input argument
            pause(1);                                                       %Pause to allow a different RNG('shuffle') state on each call
        end
        clear BCIfitResults                                                 %Free memory
    else
        if ~fileConverged(fileNrsRunningNow(newIdx))                        %If this is the first time that this file was checked (it may have been chosen multiple times, but we want the message only once)
            disp(['File nr ' num2str(fileNrsRunningNow(newIdx)) ' (' FileNames{fileNrsRunningNow(newIdx)} ') does not contain a BCIfitResults variable. This file will be skipped.']);
            fileConverged(fileNrsRunningNow(newIdx)) = 1;
            fileCounts(:,fileNrsRunningNow(newIdx)) = inf;                  %Ensure that this file is not chosen to run again
        end
        NotStartedBools(newIdx) = 1;
    end
end

%Keep calling MCMC again until we have reached converged solutions for all files in dataPath or until the number of maxAttempts has been reached for all files)  
nFinished = 0;
need2CancelBool = 0;
AllconvergedBool = 0;

while ~AllconvergedBool && any(fileCounts(2,:) < maxAttemptsPerFile)
    
    %Call MCMC again after one worker previously finished its job
    if nFinished > 0
        
        %The new index is always equal to nParallelWorkers because we deleted one of the jobs that finished   
        newIdx = nParallelWorkers;
        
        %Choose new file number (randomly one with lowest counts)
        IdxlowestCounts = find(fileCounts(1,:) == min(fileCounts(1,:)));
        newFileNr = IdxlowestCounts(randi(numel(IdxlowestCounts),1));
        fileNrsRunningNow(newIdx) = newFileNr;
        
        %Start a new job on the core that had previously finished (or whose job was cancelled)    
        StartTimes(newIdx) = tic;
        fileCounts(1,newFileNr) = fileCounts(1,newFileNr)+1;
        fileNameTmp = [dataPath filesep FileNames{fileNrsRunningNow(newIdx)}];
        listOfVarsTmp = who('-file', fileNameTmp);
        if ismember('BCIfitResults', listOfVarsTmp) 
            load(fileNameTmp,'BCIfitResults');
            if convergedAlready(BCIfitResults)
                disp(['File nr ' num2str(fileNrsRunningNow(newIdx)) ' (' FileNames{fileNrsRunningNow(newIdx)} ') has already converged. This file will be skipped.']);
                NotStartedBools(newIdx) = 1;
                fileConverged(fileNrsRunningNow(newIdx)) = 1;
                fileCounts(:,fileNrsRunningNow(newIdx)) = inf;                          %Ensure that this file is not chosen to run again
            else
                modifyBCIfitResults;                                                    %Modify some settings in BCIfitResults (see nested helper function below)
                FuturesArray(newIdx) = parfeval(@BCIfitModel,1,BCIfitResults);          %Expect 1 output argument, and use 1 input argument
                pause(1);                                                               %Pause to ensure a different RNG('shuffle') state on each call
            end
            clear BCIfitResults                                                         %Free memory
        else
            disp(['File nr ' num2str(fileNrsRunningNow(newIdx)) ' (' FileNames{fileNrsRunningNow(newIdx)} ') does not contain a BCIfitResults variable. This file will be skipped.']);
            NotStartedBools(newIdx) = 1;
            fileConverged(fileNrsRunningNow(newIdx)) = 1;
            fileCounts(:,fileNrsRunningNow(newIdx)) = inf;                              %Ensure that this file is not chosen to run again
        end
    end
    
    %Update the counter (prematurely)
    nFinished = nFinished+1;
    
    %Can we cancel any of the ongoing computations? 
    %Only cancel attempts when this file has converged already, do not cancel if it is still ongoing but when maxAttempts has been reached for this file (but not yet for all files)..    
    idx_runningConverged = find(ismember(fileNrsRunningNow,find(fileConverged)),1);
    if ~isempty(idx_runningConverged)
        need2CancelBool = idx_runningConverged;
    end
    
    %Cancel one of the ongoing futures if it has the same subjIdx as one of the converged solutions 
    if need2CancelBool
        if NotStartedBools(need2CancelBool)
            NotStartedBools(need2CancelBool) = 0;                           %no need to cancel the job if it has not even started (e.g. because 'BCIfitResults' was not found in the file) - reset the boolean now
        else
            disp(['Another future job also working on file nr ' num2str(fileNrsRunningNow(need2CancelBool)) ' (' FileNames{fileNrsRunningNow(need2CancelBool)} ') ' ...
                  'was cancelled (run time: ' datestr(toc(StartTimes(need2CancelBool))/(24*60*60),'dd HH:MM:SS') ').']);
            cancel(FuturesArray(need2CancelBool));
        end
        completedIdx = need2CancelBool;                                     %A new job will be started on this completedIdx (on the next loop within this 'while'-statement)
        need2CancelBool = 0;                        
    
    %Retrieve the MCMC results from the parallel workers as they become available       
    else
        [completedIdx,BCIfitResults] = fetchNext(FuturesArray(1:nParallelWorkers)); 
        fileCounts(2,fileNrsRunningNow(completedIdx)) = fileCounts(2,fileNrsRunningNow(completedIdx))+1;
        
        %Check whether the solution converged and whether we should save the result
        if BCIfitResults.MCMC.converge.exitflag > 0
            disp(['File nr ' num2str(fileNrsRunningNow(completedIdx)) ' (' FileNames{fileNrsRunningNow(completedIdx)} ') ' ...
                  'completed (attempt ' num2str(fileCounts(2,fileNrsRunningNow(completedIdx))) ') and converged in ' datestr(toc(StartTimes(completedIdx))/(24*60*60),'dd HH:MM:SS') '.']);
            fileConverged(fileNrsRunningNow(completedIdx)) = 1;
            fileCounts(:,fileNrsRunningNow(completedIdx)) = inf;            %Ensure that this file is not chosen to run again
            saveNowBool = 1;                                                %Save converged result
        else
            disp(['File nr ' num2str(fileNrsRunningNow(completedIdx)) ' (' FileNames{fileNrsRunningNow(completedIdx)} ') ' ...
                  'completed (attempt ' num2str(fileCounts(2,fileNrsRunningNow(completedIdx))) ') but did not converge (run time: ' datestr(toc(StartTimes(completedIdx))/(24*60*60),'dd HH:MM:SS') ').']);
            if isnan(fileMaxRs(fileNrsRunningNow(completedIdx)))
                saveNowBool = 1;                                            %Save non-converged result if no previous result for this subject has become available yet
            elseif BCIfitResults.MCMC.converge.maxR < fileMaxRs(fileNrsRunningNow(completedIdx)) 
                saveNowBool = 1;                                            %Save non-converged result if it is better than previous non-converged result in terms of maxR
            else
                saveNowBool = 0;                                            %Don't save this result
                clear BCIfitResults;                                        %Free memory
            end
            if fileCounts(2,fileNrsRunningNow(completedIdx)) >= maxAttemptsPerFile
                fileCounts(1,fileNrsRunningNow(completedIdx)) = inf;        %Ensure that this file is not chosen to run again
            end
        end
        
        %Save the result
        if saveNowBool
            fileMaxRs(fileNrsRunningNow(completedIdx)) = BCIfitResults.MCMC.converge.maxR;
            fileNameTmp = [dataPath filesep FileNames{fileNrsRunningNow(completedIdx)}];
            if appendBool
                save(fileNameTmp,'BCIfitResults','-append');
            else
                save(fileNameTmp,'BCIfitResults','-v7.3');
            end
            clear BCIfitResults;                                            %Free memory
        end
    end
    
    %Update administration by deleting the completed job id from all trackers   
    FuturesArray(completedIdx) = [];                                        %Release memory (https://uk.mathworks.com/matlabcentral/answers/424473-parfeval-memory-consumption-piling-up-clear-output-data)   
    fileNrsRunningNow(completedIdx) = []; 
    NotStartedBools(completedIdx) = []; 
    StartTimes(completedIdx) = []; 
    
    %Check whether all of the MCMC fits converged - if so we'll quit the while loop   
    AllconvergedBool = sum(fileConverged) == nFiles;
    
    %If maximum attempts has been reached for all files, but not all have converged, then break from the while loop  
    if ~AllconvergedBool && all(fileCounts(2,:) >= maxAttemptsPerFile)
        disp(['MCMC fit did converge for all files (in max = ' num2str(maxAttemptsPerFile) ' attempts per file)! Check the data and try again later using different settings...']);
        break; %break from while loop
    end
    
end %End of while loop (All MCMC converged or maxRuns reached?)

%Cancel any remaining jobs on the cores (running or queued)
cancel(FuturesArray);
delete(FuturesArray);       %Delete the objects
clear FuturesArray          %Clear the references to those objects (https://uk.mathworks.com/help/parallel-computing/parallel.job.delete.html)    

%Clear the parallel pool (to get rid of possible lingering memory leaks: https://uk.mathworks.com/matlabcentral/answers/332792-how-to-avoid-memory-leaks-when-function-inside-parfor-generates-warning)
delete(currPool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This nested function can access the workspace of the parent function
%It is memory efficient because there is no need to copy variables into input/output arguments.   
function modifyBCIfitResults

%Update settings to perform MCMC only
BCIfitResults.settings.BADS = 0; 
BCIfitResults.settings.VBMC = 0; 
BCIfitResults.settings.MCMC = 1; 
BCIfitResults.settings.Display = 0; 
BCIfitResults.settings.parallel.MCMC = 0; 
BCIfitResults.settings.nAttempts.MCMC = [1 8];  % max = 8 => thinning level = 2^7 = 128

%Continue an old non-converged MCMC run or start a new run with new initializations?  
if isfield(BCIfitResults,'MCMC')
    if isnan(fileMaxRs(fileNrsRunningNow(newIdx)))
        fileMaxRs(fileNrsRunningNow(newIdx)) = BCIfitResults.MCMC.converge.maxR;                    %Register the already existing maxR
    end
    if BCIfitResults.MCMC.converge.maxR < 1.5                                                       %If convergence is reasonably near, 
        BCIfitResults.settings.nAttempts.MCMC(2) = BCIfitResults.settings.nAttempts.MCMC(2) + 1;    %then attempt to reach convergence by sampling more.
    else                                                                                            %Otherwise,
        BCIfitResults = rmfield(BCIfitResults,'MCMC');                                              %Start a new MCMC attempt...
    end
end

end %[End of Nested Function]

end %[End of Main/Parent Function]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Function %%%
%%%%%%%%%%%%%%%%%%%%%%%

%This function cannot access the workspace of the main function, but since it does not modify its input argument no additional memory is necessary   
function bool = convergedAlready(BCIfitResults)

bool = false;
if isfield(BCIfitResults,'MCMC')
    if BCIfitResults.MCMC.converge.exitflag > 0
        bool = true;
    end
end

end %[EoF]
