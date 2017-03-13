% RR3 group time-frequency plots (eeglab)
% Created: 06/25/2015 By: Evan Layher
% Revised: 11/14/2015 By: Evan Layher % Prompt which group files to get
% Revised: 12/15/2015 By: Evan Layher (input single or grouped graphs)
% Adapted from Frank's code: /nfs/to-eeg/code/frank_scripts

% Creates time frequency plots using tftopo()
% Working Memory (wm) and Temporal Sequence (ts)

clear all
clc

addpath('/nfs/pkg64/eeglab/eeglab13_0_1b/');  % eeglab directory
addpath('/nfs/to-eeg/code/');                 % code directory

% WORKING MEMORY VALUES
wmDataPath  = '/nfs/to-eeg/wm_data/group/';
wmVertLines = [-9500, -7500, -7000, -5000, -4500, -2500, -2000, 0];
wmFreqRange = [3 40];   % input freqs (Hz)
wmTimeRange = [0 4000]; % input times (ms)


groupaveragedfilestoplot = {'wm_baseline_sr6_freqs_3_40_27CN_21SZmed_8SZunmed_2017_01_17_fh.mat',...
                            'wm_no_baseline_sr6_freqs_3_40_27CN_21SZmed_8SZunmed_2017_01_17_fh.mat'};
groupContrasts = {'combined.cnAvg' 'combined.szAvg' 'combined.cnVsz' 'combined.cnVsznomeds' 'combined.cnVszsznomeds' 'combined.szsznomedsAvg'}; % List of group structures
noBaselineName = '_no_baseline_'; % Part of filename that specifies no baseline data


% Graph options
cRangeContrast = [-2 2];  % Blue to Red color scale for contrast images
cRangeBase     = [-2 2];  % Blue to Red color scale (baseline subtraction)
cRangeNoBase   = [20 60]; % Blue to Red color scale (no baseline subtraction)
addVertLines   = true; % true or false (Vertical lines in graph)
addColorBar    = true; % true or false (Colorbar on graph)
outFileExtension = '.jpg';

defaultPlotWidth  = 1.86; % Inches (1.86x1.40 fits 9 plots in powerpoint)
defaultPlotHeight = 1.40; % Inches (1.86x1.40 fits 9 plots in powerpoint)
defaultColorBarSettings = [0.9, 0.15, 0.02, 0.7]; % Percentages [Left->Right, Bottom->Top, width, height]
defaultGraphPosition = [0.08, 0.12, 0.8, 0.77]; % Percentages [Left->Right, Bottom->Top, width, height]
defaultFontSize = 6; 

%%%%%-----USER INPUT OPTIONS-----%%%%%
fprintf('CREATE RR3 AVERAGE GROUP TIME-FREQUENCY PLOTS\n\n')
% INPUT TASK TYPE(s) (Working memory and/or temporal sequence)
validTask = true; % Continues while loop until false
while validTask
    fprintf('WHICH TASK(s)?\n')
    taskPrompt = 'input [1]: wm AND ts time-frequency plots\ninput [2]: ts time-frequency plots\ninput [3]: wm time-frequency plots:';
    taskInput = input(taskPrompt);
    if taskInput == 1
        task = {'ts' 'wm'};
        validTask = false;
    elseif taskInput == 2
        task = {'ts'};
        validTask = false;
    elseif taskInput == 3
        task = {'wm'};
        validTask = false;
    else
        clc
        taskInput = num2str(taskInput);
        fprintf('INVALID INPUT: %s\n', taskInput)
    end
end

% INPUT GRAPH TYPE(s) (single and/or grouped electrode graphs)
validGraph = true; % Continues while loop until false
createSingleElectrodeGraphs = false; % true or false create single electrode graphs
createGroupElectrodeGraphs  = false; % true or false create grouped electrode graphs
while validGraph
    fprintf('\nWHICH ELECTRODE PLOT TYPE(s)?\n')
    graphPrompt = 'input [1]: grouped AND single electrode time-frequency plots\ninput [2]: grouped electrode time-frequency plots\ninput [3]: single electrode time-frequency plots:';
    graphInput = input(graphPrompt);
    if graphInput == 1
        createSingleElectrodeGraphs = true; % true or false create single electrode graphs
        createGroupElectrodeGraphs  = true; % true or false create grouped electrode graphs
        validGraph = false;
    elseif graphInput == 2
        createGroupElectrodeGraphs  = true; % true or false create grouped electrode graphs
        validGraph = false;
    elseif graphInput == 3
        createSingleElectrodeGraphs = true; % true or false create single electrode graphs
        validGraph = false;
    else
        clc
        graphInput = num2str(graphInput);
        fprintf('INVALID INPUT: %s\n', graphInput)
    end
end

% INPUT OUTPUT PICTURE DIMENSIONS
validDimensions = true;
while validDimensions
    fprintf('\nWHAT DIMENSIONS FOR OUTPUT PLOTS?\n')
    fprintf('Default time-frequency plot dimensions in inches: %.2fin x %.2fin\n', defaultPlotWidth, defaultPlotHeight)
    dimensionPrompt = 'input [0] to use default inch dimensions\ninput [width height] for custom dimensions:';
    dimensionInput = input(dimensionPrompt);
    
    if dimensionInput == 0
        plotWidth  = defaultPlotWidth;
        plotHeight = defaultPlotHeight;
        validDimensions = false;
    elseif isfloat(dimensionInput) && length(dimensionInput) == 2 && dimensionInput(1) > 0 && dimensionInput(1) < 100 && dimensionInput(2) > 0 && dimensionInput(2) < 100
        plotWidth  = dimensionInput(1);
        plotHeight = dimensionInput(2);
        validDimensions = false;
    else
        clc
        dimensionInput = num2str(dimensionInput);
        fprintf('INVALID INPUT: %s\n', dimensionInput)
    end
end

%%%%%-----START SCRIPT-----%%%%%
tStart = tic; % Script start time
close all     % Close open matlab windows
eeglab        % load eeglab first

%%% Group electrodes (Variable names will become titles in graphs)
%PFC left
Left_PFC = [1 4 6 7 8 15 16 17];

%PFC mid
Mid_PFC = [2 9 10 11 18 19 20];

%PFC right
Right_PFC = [3 5 12 13 14 21 22 23];

%Central left
Left_Central = [24 25 26 34 35 36];

%Central mid
Mid_Central = [27 28 29 37 38 39];

%Central right
Right_Central = [30 31 32 40 41 42];

%Posterior left
Left_Post = [44 45 46 53 54];

%Posterior mid
Mid_Post = [47 48 49 55 56 57 61 62 63];

%Posterior right
Right_Post = [50 51 52 58 59];

% Input region variables as 'strings'
groupRegions = {'Left_PFC', 'Mid_PFC', 'Right_PFC', 'Left_Central', 'Mid_Central', 'Right_Central', 'Left_Post', 'Mid_Post', 'Right_Post'};

for iTask = 1:length(task) % wm and/or ts
    eegTask = task{iTask};
    if strcmp(eegTask, 'ts')
        dataPath  = tsDataPath;
        vertLines = tsVertLines;
        timeRange = tsTimeRange;
        freqRange = tsFreqRange;
    elseif strcmp(eegTask, 'wm')
        dataPath  = wmDataPath;
        vertLines = wmVertLines;
        timeRange = wmTimeRange;
        freqRange = wmFreqRange;
    else
        fprintf('ERROR: task must be "wm" and/or "ts"\nTask currently equals: "%s"\n', eegTask)
        continue
    end
    
%     dataFiles = dir(dataPath);
%     askCount  = 0; % Build structure of potential group files
%     
%     for jAskFile = 1:length(dataFiles) % Loop through each file found
%         checkDataFile = [dataPath, dataFiles(jAskFile).name];
%         if strcmp(checkDataFile(end-3:end), '.mat') % Confirm matfile
%             askCount = askCount + 1;
%             askFile(askCount).name = checkDataFile;
%             fprintf('[%d] %s\n', askCount, askFile(askCount).name) % Display for user input
%         else
%             fprintf('Not a ".mat" file: %s\n', checkDataFile)
%         end
%     end
%     
%     if askCount == 0
%         error('NO GROUP ".mat" FILES FOUND: %s\n', dataPath)
%     elseif askCount < 3 % 1 group type found (baseline and/or no baseline)
%         finalFiles = askFile;
%     else % Prompt user to specify which group '.mat' files to use
%         validIndex = true; % Continues while loop until false
%         while validIndex   % Must input value from 1 to length(askFile)
%             indexPrompt = 'Input numbers of desired group ".mat" files to get time-frequency plots from (e.g. [1 6]):';
%             indices     = input(indexPrompt); % Subject ID 
%             if isempty(indices) % Must specify indices
%                 fprintf('NO INPUTS SPECIFIED\n')
%             else % Make sure valid input
%                 validIndex = false;
%                 finalCount = 0;
%                 for jIndex = indices
%                     finalCount = finalCount + 1;
%                     if jIndex >= 1 && jIndex <= askCount % Index must be within range
%                         finalFiles(finalCount).name = askFile(jIndex).name;
%                     else
%                         fprintf('INVALID INPUT: %d\n', jIndex)
%                         validIndex = true;
%                         break    
%                     end
%                 end
%             end
%         end
%     end
%     
%     for jDisplayFile = 1:length(finalFiles) % Display files to plot
%         dataDisplayFile = finalFiles(jDisplayFile).name;
%         fprintf('CREATING PLOTS: %s\n', dataDisplayFile)
%     end
    
    for jFile = 1:length(groupaveragedfilestoplot)
        dataFile = [dataPath, groupaveragedfilestoplot{jFile}];
        load(dataFile);
        if ~isempty(strfind(dataFile, noBaselineName))
            analysisType = '_no_baseline_';
        else
            analysisType = '_baseline_';
        end
        
        
        % Output directories
        outPath      = [dataPath, 'tf_plots/'];
        outPathElec  = [outPath, 'single_elecs/'];
        outPathGroup = [outPath, 'group_elecs/'];
        
        if ~exist(outPath, 'dir')
            mkdir(outPath)
        end
        
        % get frequency range
        if isempty(freqRange) % Use frequency information in files
            freqs = combined.freqs;
            freqPlot = 1:length(combined.freqs);
        else % Use specified frequency range
            freqFind = find(combined.freqs == freqRange(1) | combined.freqs == freqRange(2));
            if length(freqFind) == 2
                freqPlot = freqFind(1):freqFind(2);
                freqs = combined.freqs(freqFind(1):freqFind(2));
            else
                error('Invalid frequency range: %.2f - %.2f Hz FILE: %s\n', freqRange(1), freqRange(2), dataFile)
            end
        end
        
        % get time range
        if isempty(timeRange)
            times = combined.times; % Use time information in files
            timePlot = 1:length(combined.times);
        else % Use specified time range
            timeFind = find(combined.times == timeRange(1) | combined.times == timeRange(2));
            if length(timeFind) == 2
                timePlot = timeFind(1):timeFind(2);
                times = combined.times(timeFind(1):timeFind(2));
            else
                error('Invalid time range: %d - %d ms Subject: %s\n', timeRange(1), timeRange(2), dataFile)
            end
        end
        
        % Omit vertical lines if outside of time range
        if ~isempty(vertLines)
            if vertLines(1) < times(1) || vertLines(end) > times(end)
                vertLines = [];
                fprintf('\n\nEXCLUDING VERTICAL LINES: TIME RANGE TOO SHORT\n\n')
            end
        end
        
        % Include frequency and time ranges in filename
        freqTimeOutputName = [num2str(freqs(1)), '_', num2str(freqs(end)), 'Hz_', num2str(times(1)), '_', num2str(times(end)), 'ms', outFileExtension];
        
        % Group contrasts
        for kContrast = 1:length(groupContrasts)
            contrastType = groupContrasts{kContrast};
            contrastName = strrep(contrastType, 'combined.', '');
            % trial types
            for mTrial = 1:length(combined.conditions)
                trialName = combined.conditions{mTrial};
                
                if ~isempty(strfind(trialName,'_v_')) || isempty(strfind(contrastName),'V') % Order vs Item Hit || cnVsz 
                    cRange = cRangeContrast;
                elseif strcmp(analysisType, '_no_baseline_')
                    cRange = cRangeNoBase;
                else
                    cRange = cRangeBase;
                end
                
                % Individual electrodes
                if createSingleElectrodeGraphs % true/false
                    if ~exist(outPathElec, 'dir')
                        mkdir(outPathElec)
                    end
                    
                    for nElec = 1:length(combined.chanLabels) % Loop through electrodes
                        graphTitle = [contrastName, '_', trialName, '_', combined.chanLabels{nElec}];
                        close all % close all figures
                        if addVertLines % true/false
                            graphPlot = tftopo(eval([contrastType, '.conditions{', num2str(mTrial), '}.chan{', num2str(nElec), '}(freqPlot,timePlot)']), times, freqs, 'vert', vertLines);
                        else
                            graphPlot = tftopo(eval([contrastType, '.conditions{', num2str(mTrial), '}.chan{', num2str(nElec), '}(freqPlot,timePlot)']), times, freqs);
                        end
                        set(gca, 'FontSize', defaultFontSize);      % axis number font size
                        set(gca, 'Position', defaultGraphPosition); % graph position
                        xlabel('ms', 'FontSize', defaultFontSize);  % x-axis label
                        ylabel('Hz', 'FontSize', defaultFontSize);  % y-axis label
                        title(strrep(graphTitle, '_', ' '), 'FontSize', defaultFontSize);
                        caxis(cRange); % Power Range
                        set(gcf, 'PaperUnits', 'inches'); % Input units of inches
                        set(gcf, 'PaperPosition', [0, 0, plotWidth, plotHeight]); % Graphic Size
                        graphOutputName = [graphTitle, analysisType, freqTimeOutputName];
                        graphOutputFile = [outPathElec, graphOutputName];
                        if addColorBar % true/false
                            set(colorbar, 'Position', defaultColorBarSettings, 'FontSize', defaultFontSize);
                        end
                        saveas(gcf, graphOutputFile);
                        fprintf('CREATED: %s\n', graphOutputFile)
                    end % nElec
                end % if createSingleElectrodeGraphs
                
                % Groups of electrodes
                if createGroupElectrodeGraphs % true/false
                    if ~exist(outPathGroup, 'dir')
                        mkdir(outPathGroup)
                    end
                    
                    groupAverage = cell(size(groupRegions));
                    
                    figure;                                     % Added 02/16/2017 By FH
                    set(gcf, 'Position',get(0,'Screensize'));   % Added 02/16/2017 By FH
                    for nRegion = 1:length(groupRegions)
                        % groupRegions are strings for naming files
                        groupReg = eval(groupRegions{nRegion}); % (e.g.: groupReg = Left_PFC  (Left_PFC = [1 4 6 7 8 15 16 17])
                        elecArray = nan(length(freqs), length(times), length(groupReg));
                        
                        for nElec = 1:length(groupReg) % Loop through each electrode in region
                            elecArray(:, :, nElec) = eval([contrastType, '.conditions{', num2str(mTrial), '}.chan{', num2str(groupReg(nElec)), '}(freqPlot,timePlot)']);
                        end % (e.g.: elecArray(:,:,nElec) = combined.cnAvg.trial{1}.chan{50}(freqPlot,timePlot)
                        
                        groupAverage{nRegion} = mean(elecArray, 3); % Average region electrodes
                        graphTitle = [contrastName, '_', trialName, '_', groupRegions{nRegion}];
%                       commendted out 02/16/2015 By FH  close all
                        
                        if addVertLines % true/false
                            subplot(3,3,nRegion);                                               % Added 02/16/2017 By FH
                            tftopo(groupAverage{nRegion}, times, freqs, 'vert', vertLines);     % Added 02/16/2017 By FH
%                           commented out 02/16/2017 By FH  graphPlot = tftopo(groupAverage{nRegion}, times, freqs, 'vert', vertLines);
                        else
                            subplot(3,3,nRegion);                                               % Added 02/16/2017 By FH
                            tftopo(groupAverage{nRegion}, times, freqs);                        % Added 02/16/2017 By FH
%                           commented out 02/16/2017 By FH  graphPlot = tftopo(groupAverage{nRegion}, times, freqs);
                        end
                        set(gca, 'FontSize', defaultFontSize);      % axis number font size
%                       commented out 02/16/2017 By FH  set(gca, 'Position', defaultGraphPosition); % graph position
                        xlabel('ms', 'FontSize', defaultFontSize);  % x-axis label
                        ylabel('Hz', 'FontSize', defaultFontSize);  % y-axis label
                        title(strrep(graphTitle, '_', ' '), 'FontSize', defaultFontSize);
                        caxis(cRange); % Power Range
%                       commented out 02/16/2017 By FH  set(gcf, 'PaperUnits', 'inches'); % Input units of inches
%                       commented out 02/16/2017 By FH  set(gcf, 'PaperPosition', [0, 0, plotWidth, plotHeight]); % Graphic Size
%                       commented out 02/16/2017 By FH  graphOutputName = [graphTitle, '_CN', analysisType, freqTimeOutputName];
%                       commented out 02/16/2017 By FH  graphOutputFile = [outPathGroup, graphOutputName];
                        if addColorBar % true/false
%                           commented out 02/16/2017 By FH  set(colorbar, 'Position', defaultColorBarSettings, 'FontSize', defaultFontSize);
                            colorbar;       % Added 02/16/2017 By FH
                        end
%                       commented out 02/16/2017 By FH  saveas(gcf, graphOutputFile);
%                       commented out 02/16/2017 By FH  fprintf('CREATED: %s\n', graphOutputFile)
                    end % mRegion
                    graphOutputName = [contrastName, '_', trialName, analysisType, freqTimeOutputName]; % Added 02/16/2017 By FH
                    graphOutputFile = [outPathGroup, graphOutputName];                                  % Added 02/16/2017 By FH
                    saveas(gcf, graphOutputFile);                                                       % Added 02/16/2017 By FH
                    fprintf('CREATED: %s\n', graphOutputFile);                                          % Added 02/16/2017 By FH
                    close all;                                                                          % Added 02/16/2017 By FH
                end % if createGroupElectrodeGraphs
            end % mTrial
        end % kContrast
    end % jFile
end % iTask

fprintf('PROCESS TIME: ')
hourToc(tStart); % Display subject process time