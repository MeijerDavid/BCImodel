function fig_handles = BCIdisplay(BCIfitResults)
%Quick and dirty plotting function. Create somes figures to compare 
%response distributions to BCI-predicted likelihood functions. 
%This function can also be called directly after fitting has finished...

fig_handles = [];

%If this function was called only to display (BADS, VBMC and MCMC were not executed) - then display everything!
if ~BCIfitResults.settings.BADS && ~BCIfitResults.settings.VBMC && ~BCIfitResults.settings.MCMC
    BCIfitResults.settings.BADS = 1; BCIfitResults.settings.VBMC = 1; BCIfitResults.settings.MCMC = 1;
end

%Find all location combinations
sAsVconds = BCIfitResults.settings.sAsVconds;
uniqA = unique(sAsVconds(~isnan(sAsVconds(:,1)),1))';
nAlocsAresp = numel(uniqA);
if any(isnan(sAsVconds(:,1))); uniqA = [uniqA NaN]; end
uniqV = unique(sAsVconds(~isnan(sAsVconds(:,2)),2))';
nVlocsVresp = numel(uniqV);
if any(isnan(sAsVconds(:,2))); uniqV = [uniqV NaN]; end
nAlocs = length(uniqA); nVlocs = length(uniqV);

%Find number of conditions and number of responses per condition
[nLocsComb,nConds] = size(BCIfitResults.BADS.goodnessOfFit.predictedResps);
[nAresps,nVresps,nCSJresps] = cellfun(@extractCondInfo,BCIfitResults.BADS.goodnessOfFit.predictedResps);

%Discrete responses?
if all(~isnan(BCIfitResults.settings.RespLocs)) 
    discreteBool = 1;
else
    discreteBool = 0;
end

%Response range
RespRange = BCIfitResults.settings.RespRange;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Report quantitative fitting measures in command window (BADS or MCMC) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(BCIfitResults,'MCMC')
    GoFstruct = BCIfitResults.MCMC.goodnessOfFit;
else
    GoFstruct = BCIfitResults.BADS.goodnessOfFit;
end

if BCIfitResults.settings.BADS || BCIfitResults.settings.MCMC
    disp(' '); disp('Overall quantitative fitting measures: ...')
    printQuantFittingMeasures(GoFstruct,1,discreteBool,sum(nAresps,'all'),sum(nVresps,'all'),sum(nCSJresps,'all'),1);
    %Report condition-specific quantitative fitting measures in command window
    for c=1:nConds
        if nConds > 1
            disp(' '); disp(['Condition #' num2str(c) ': quantitative fitting measures: ...']);
            printQuantFittingMeasures(GoFstruct,c+1,discreteBool,sum(nAresps(:,c)),sum(nVresps(:,c)),sum(nCSJresps(:,c)),0);
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the response distributions (BADS) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.BADS
    for c=1:nConds

        %Prepare to put R2 fit measures in figure names if more than one condition   
        if nConds > 1
            if sum(nAresps(:,c)) > 0
                A_str =   ['R2_N = ' num2str(GoFstruct.R2{2,c+1}.respA_only.BCIvsNull_Nagelkerke) '. ' ...
                           'R2_S = ' num2str(GoFstruct.R2{2,c+1}.respA_only.BCIvsNull_Stochastic) '. '];
                if discreteBool; A_str = [A_str 'AbsGoF = ' num2str(GoFstruct.absGoF{2,c+1}.respA_only)]; end      
            end
            if sum(nVresps(:,c)) > 0
                V_str =   ['R2_N = ' num2str(GoFstruct.R2{2,c+1}.respV_only.BCIvsNull_Nagelkerke) '. ' ...
                           'R2_S = ' num2str(GoFstruct.R2{2,c+1}.respV_only.BCIvsNull_Stochastic) '. '];
                if discreteBool; V_str = [V_str 'AbsGoF = ' num2str(GoFstruct.absGoF{2,c+1}.respV_only)]; end       
            end
            if sum(nCSJresps(:,c)) > 0
                CSJ_str = ['R2_N = ' num2str(GoFstruct.R2{2,c+1}.respCSJ_only.BCIvsNull_Nagelkerke) '. ' ...
                           'R2_S = ' num2str(GoFstruct.R2{2,c+1}.respCSJ_only.BCIvsNull_Stochastic) '. ' ...
                           'AbsGoF = ' num2str(GoFstruct.absGoF{2,c+1}.respCSJ_only)];         
            end
        else
            A_str = []; V_str = []; CSJ_str = [];
        end
        
        %Find the amount of motor noise (fitted or fixed)
        if ~discreteBool
            if any(strcmp(BCIfitResults.settings.ParamNames2Fit,'sigmaMotor'))
                paramidx = find(strcmp(BCIfitResults.settings.ParamNames2Fit,'sigmaMotor'));
                paramidxC = paramidx(ismember(paramidx,BCIfitResults.settings.ParamsPerCond{c}));
                sigmaMotor = BCIfitResults.BADS.fittedParams(paramidxC(1));
            else
                sigmaMotor = BCIfitResults.settings.sigmaMotor;
            end
        end
        
        %Plot Auditory Response Distributions
        if sum(nAresps(:,c)) > 0 

            %Open a new figure
            tempTitle = ['AUD_ Auditory Responses: condition #' num2str(c) '. BCI = red, Data = blue. ' A_str]; 
            handleTmp = figure('units','normalized','outerposition',[0 0 1 1],'NumberTitle', 'off', 'Name', tempTitle);
            fig_handles = [fig_handles handleTmp];
            
            for i=1:nLocsComb

                sA = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sA;
                if isnan(sA)
                    iA = nAlocs;
                else
                    iA = find(uniqA == sA);
                end

                sV = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sV;
                if isnan(sV)
                    iV = nVlocs;
                else
                    iV = find(uniqV == sV);
                end

                %Don't plot auditory responses if sA = NaN
                if ~(any(isnan(uniqA)) && (iA == nAlocs))

                    subPlotIdx = (iV-1)*nAlocsAresp+iA;
                    subplot(nVlocs,nAlocsAresp,subPlotIdx); box on; hold on;

                    %Plot BCI prediction
                    xBCI = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.likeFunA_grid;
                    yBCI = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.likeFunA;

                    plot([sV sV],[0 max(yBCI)],'y--','Linewidth',2);
                    plot([sA sA],[0 max(yBCI)],'g:','Linewidth',2);
                    plot(xBCI,yBCI,'r');

                    %Plot normalized histogram of responses
                    Aresps = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.respA';
                    Aresps = Aresps(~isnan(Aresps));
                    if discreteBool
                        %edges = [-inf xBCI(1:(end-1))+diff(xBCI)/2 inf];   
                        edges = (0:length(xBCI))+0.5;                       %Responses were already discretized    
                        yResp = histcounts(Aresps,edges,'Normalization','probability');
                        plot(xBCI,yResp,'b');
                    else
                        nBins = ceil(numel(xBCI)/3);
                        nAtemp = numel(Aresps);
                        [yResp,xResp] = BCIcompLikeFun(sort(Aresps),ones(1,nAtemp)/nAtemp,nBins,RespRange,2,sigmaMotor);
                        plot(xResp,yResp,'b');
                    end
                    
                    title(['sA = ' num2str(sA) '°, sV = ' num2str(sV) '°']);
                    if (iV==nVlocs); xlabel('Response location (°)'); end
                    if (iA==1)
                        if discreteBool
                            ylabel('Likelihood (prob)');
                        else
                            ylabel('Likelihood (pdf)');
                        end
                    end
                    xlim(RespRange);
                    if discreteBool
                        xticks(BCIfitResults.settings.RespLocs);
                    end
                end
            end
        end 

        %Plot Visual Response Distributions
        if sum(nVresps(:,c)) > 0 

            %Open a new figure
            tempTitle = ['VIS_ Visual Responses: condition #' num2str(c) '. BCI = red, Data = blue. ' V_str]; 
            handleTmp = figure('units','normalized','outerposition',[0 0 1 1],'NumberTitle', 'off', 'Name', tempTitle);
            fig_handles = [fig_handles handleTmp];
            
            for i=1:nLocsComb

                sA = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sA;
                if isnan(sA)
                    iA = nAlocs;
                else
                    iA = find(uniqA == sA);
                end

                sV = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sV;
                if isnan(sV)
                    iV = nVlocs;
                else
                    iV = find(uniqV == sV);
                end

                %Don't plot visual responses if sV = NaN
                if ~(any(isnan(uniqV)) && (iV == nVlocs))

                    subPlotIdx = (iV-1)*nAlocs+iA;
                    subplot(nVlocsVresp,nAlocs,subPlotIdx); box on; hold on;

                    %Plot BCI prediction
                    xBCI = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.likeFunV_grid;
                    yBCI = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.likeFunV;

                    plot([sA sA],[0 max(yBCI)],'g--','Linewidth',2);
                    plot([sV sV],[0 max(yBCI)],'y:','Linewidth',2);
                    plot(xBCI,yBCI,'r');

                    %Plot normalized histogram of responses
                    Vresps = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.respV';
                    Vresps = Vresps(~isnan(Vresps));
                    if discreteBool
                        %edges = [-inf xBCI(1:(end-1))+diff(xBCI)/2 inf];
                        edges = (0:length(xBCI))+0.5;                       %Responses were already discretized
                        yResp = histcounts(Vresps,edges,'Normalization','probability');
                        plot(xBCI,yResp,'b');
                    else
                        nBins = ceil(numel(xBCI)/3);
                        nVtemp = numel(Vresps);
                        [yResp,xResp] = BCIcompLikeFun(sort(Vresps),ones(1,nVtemp)/nVtemp,nBins,RespRange,2,sigmaMotor);
                        plot(xResp,yResp,'b');
                    end

                    title(['sA = ' num2str(sA) '°, sV = ' num2str(sV) '°']);
                    if (iV==nVlocsVresp); xlabel('Response location (°)'); end
                    if (iA==1)
                        if discreteBool
                            ylabel('Likelihood (prob)');
                        else
                            ylabel('Likelihood (pdf)');
                        end
                    end
                    xlim(RespRange);
                    if discreteBool
                        xticks(BCIfitResults.settings.RespLocs);
                    end
                end
            end
        end 

        %Plot Common-Source Judgement Response Distributions
        if sum(nCSJresps(:,c)) > 0 

            %Open a new figure
            tempTitle = ['CSJ_ Common Source Judgments: condition #' num2str(c) '. BCI = red, Data = blue. ' CSJ_str]; 
            handleTmp = figure('units','normalized','outerposition',[0 0 1 1],'NumberTitle', 'off', 'Name', tempTitle);
            fig_handles = [fig_handles handleTmp];

            for i=1:nLocsComb

                sA = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sA;
                if isnan(sA)
                    iA = nAlocs;
                else
                    iA = find(uniqA == sA);
                end

                sV = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.sV;
                if isnan(sV)
                    iV = nVlocs;
                else
                    iV = find(uniqV == sV);
                end

                %Don't plot common-source judgments if sA = NaN or if sV = NaN   
                if ~((any(isnan(uniqA)) && (iA == nAlocs)) || (any(isnan(uniqV)) && (iV == nVlocs)))

                    subPlotIdx = (iV-1)*nAlocsAresp+iA;
                    subplot(nVlocsVresp,nAlocsAresp,subPlotIdx); box on; hold on;

                    %Plot BCI prediction
                    xBCI = [1 2];
                    yBCI = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.likeFunCSJ;
                    plot(xBCI,yBCI,'r');

                    %Plot normalized histogram of responses
                    CSJresps = BCIfitResults.BADS.goodnessOfFit.predictedResps{i,c}.respCSJ';
                    CSJresps = CSJresps(~isnan(CSJresps));
                    yResp = [NaN NaN];
                    yResp(1) = sum(CSJresps == 1)/numel(CSJresps);                      %A "1" means common-source
                    yResp(2) = sum(CSJresps ~= 1 & ~isnan(CSJresps))/numel(CSJresps);   %Any other (e.g. "0" or "2") but not NaN means independent sources!
                    plot(xBCI,yResp,'b');

                    title(['sA = ' num2str(sA) '°, sV = ' num2str(sV) '°']);
                    if (iV==nVlocsVresp); xlabel('Number of sources'); end 
                    if (iA==1); ylabel('Likelihood (prob)'); end
                    xlim([0.5 2.5]);
                    ylim([0 1]);
                    xticks([1 2]);
                end
            end
        end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot Posterior Parameter Distributions (VBMC & MCMC) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if BCIfitResults.settings.VBMC && isfield(BCIfitResults,'VBMC')
    %Add relevant 'VBMC' folders to the path (just in case)
    me = mfilename;                                             
    pathstr = fileparts(fileparts(which(me)));                  
    addpath(genpath([pathstr filesep 'vbmc-master']));  
    
    Xs = BCIconvertParams2Fit(vbmc_rnd(BCIfitResults.VBMC.vp,1e6,1,1),BCIfitResults.settings,'conv2real');              %Generate a million samples from the variational posterior 
    h_fig_VBMC = cornerplot(Xs,BCIfitResults.settings.ParamNames2Fit,BCIfitResults.BADS.fittedParams);                  %Plot VBMC results in cornerplot
    set(h_fig_VBMC, 'NumberTitle', 'off', 'Name', 'VBMC - posterior parameter distributions. Black = MAPs');
    set(h_fig_VBMC,'units','normalized','outerposition',[0 0 1 1]);
    fig_handles = [fig_handles h_fig_VBMC];
end

if BCIfitResults.settings.MCMC && isfield(BCIfitResults,'MCMC')
    %Add relevant 'MCMC' subfolders to the path (just in case)
    me = mfilename;                                                 
    pathstr = fileparts(fileparts(which(me)));                      
    addpath(genpath([pathstr filesep 'eissample-master']));   
    addpath(genpath([pathstr filesep 'vbmc-master']));                      %Required by cornerplot
    
    bounds = quantile(BCIfitResults.MCMC.mcmcParams,[0.05 0.95],1);         %posterior quantiles (0.5 = median). This returns a matrix with size [2 x nParams]
    samples = 1:size(BCIfitResults.MCMC.mcmcParams,1);                      %change this if you want to look at a subsection of all samples
    h_fig_MCMC = cornerplot_DM(BCIfitResults.MCMC.mcmcParams(samples,:),BCIfitResults.settings.ParamNames2Fit,BCIfitResults.BADS.fittedParams,bounds,[],BCIfitResults.MCMC.prob.mcmcLogPosterior(samples,:));
    set(h_fig_MCMC, 'NumberTitle', 'off', 'Name', 'MCMC - posterior parameter distributions. Black = MAPs. Blue = Parameter Distributions. Red = Marginalized Probability Distributions');
    set(h_fig_MCMC,'units','normalized','outerposition',[0 0 1 1]);
    fig_handles = [fig_handles h_fig_MCMC];
end

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Note that it would be very easy to extend this function and obtain more info    
function [nAresps,nVresps,nCSJresps] = extractCondInfo(predictedRespsStruct)
    sA = predictedRespsStruct.sA;
    sV = predictedRespsStruct.sV;
    respA = predictedRespsStruct.respA;
    respV = predictedRespsStruct.respV;
    respCSJ = predictedRespsStruct.respCSJ;
    likeFunA = predictedRespsStruct.likeFunA;
    likeFunV = predictedRespsStruct.likeFunV;
    likeFunCSJ = predictedRespsStruct.likeFunCSJ;
    likeFunA_grid = predictedRespsStruct.likeFunA_grid;
    likeFunV_grid = predictedRespsStruct.likeFunV_grid;
    
    nAresps = sum(~isnan(respA));
    nVresps = sum(~isnan(respV));
    nCSJresps = sum(~isnan(respCSJ));
end

%-------------------------------------------------------------------------

function printQuantFittingMeasures(GoFstruct,iCond,discreteBool,nAresps,nVresps,nCSJresps,specificsBool)

    disp(['R2_total_Nagelkerke: ' num2str(GoFstruct.R2{2,iCond}.total.BCIvsNull_Nagelkerke)]);
    disp(['R2_total_Stochastic: ' num2str(GoFstruct.R2{2,iCond}.total.BCIvsNull_Stochastic)]);
    if discreteBool 
        if (nAresps > 0) && (nVresps > 0)
            disp(['AbsGoF_Aud&Vis: ' num2str(GoFstruct.absGoF{2,iCond}.respA_and_respV)]); 
        elseif (nAresps > 0) && ~specificsBool
            disp(['AbsGoF_AudOnly: ' num2str(GoFstruct.absGoF{2,iCond}.respA_only)]);
        elseif (nVresps > 0) && ~specificsBool
            disp(['AbsGoF_VisOnly: ' num2str(GoFstruct.absGoF{2,iCond}.respV_only)]);
        end
    end
    if specificsBool
        if nAresps > 0
            disp(['R2_AudOnly_Nagelkerke: ' num2str(GoFstruct.R2{2,iCond}.respA_only.BCIvsNull_Nagelkerke)]);
            disp(['R2_AudOnly_Stochastic: ' num2str(GoFstruct.R2{2,iCond}.respA_only.BCIvsNull_Stochastic)]);
            if discreteBool; disp(['AbsGoF_AudOnly: ' num2str(GoFstruct.absGoF{2,iCond}.respA_only)]); end
        end
        if nVresps > 0
            disp(['R2_VisOnly_Nagelkerke: ' num2str(GoFstruct.R2{2,iCond}.respV_only.BCIvsNull_Nagelkerke)]);
            disp(['R2_VisOnly_Stochastic: ' num2str(GoFstruct.R2{2,iCond}.respV_only.BCIvsNull_Stochastic)]);
            if discreteBool; disp(['AbsGoF_VisOnly: ' num2str(GoFstruct.absGoF{2,iCond}.respV_only)]); end
        end
        if nCSJresps > 0
            disp(['R2_CSJonly_Nagelkerke: ' num2str(GoFstruct.R2{2,iCond}.respCSJ_only.BCIvsNull_Nagelkerke)]);
            disp(['R2_CSJonly_Stochastic: ' num2str(GoFstruct.R2{2,iCond}.respCSJ_only.BCIvsNull_Stochastic)]);
        end
    end
    if nCSJresps > 0
        disp(['AbsGoF_CSJonly: ' num2str(GoFstruct.absGoF{2,iCond}.respCSJ_only)]);
    end
    
end %[EoF]
