% RR3 WM ERSP (eeglab)
% Created: 04/22/2015 By: Evan Layher
% Revised: 08/13/2015 By: Evan Layher
% Revised: 05/09/2016 By: Frank Hsieh
% Adapted from Frank's code: /nfs/to-eeg/code/frank_scripts

% Creates ERSP, ITC and single trial time/frequency decompositions
% ERSP is baseline subtracted from the mean baseline of individuals trials
% Creates baseline and non-baseline subtracted values in single mat file

clear all
clc

addpath('/nfs/pkg64/eeglab/eeglab13_0_1b/');  % eeglab directory
addpath('/nfs/to-eeg/code/');                 % code directory

dataPath = '/nfs/to-eeg/wm_data/';
dataFileEnding = '_7_wm_final.set';

tfCycles           = 6;      % (morlet cycles) (use 4 5 or 6)
freqRange          = [3 40]; % (Hz) [lowfreq highfreq]
tfBaseline         = [-10500 -10000]; % baseline (500ms prior to first stimulus)
includedElectrodes = 1:64;   % 64 channel neuroscan cap

lowResTimeArray    = -14500:10:5000; % (ms) output times (start:increment:end)
highResTimeArray   = -14500:2:5000;  % (ms) output times (start:increment:end)

lowResFreqInterval  = length(freqRange(1):freqRange(2)); % length([freqRange(1):freqRange(2)]) = 1 Hz interval.
HighResFreqInterval = 4 * length(freqRange(1):freqRange(2)) - 3; % Multiply by 4 (and subtract 3) to get 0.25 Hz

validSubId = true; % Continues while loop until false
while validSubId   % Must input value between 100 and 400
    subPrompt = 'Input IDs in brackets [101 102 103] or [0] for ALL IDs:';
    subjects = input(subPrompt);  % RR3 subject ID(s)
    validSubId = false;  % exit loop
    for iSubInput = subjects
        if iSubInput == 0
            subjects = []; % reset subjects array
            allDirs = dir(dataPath);
            for jDir = 1:length(allDirs)
                if iSubInput > 100 || iSubInput < 400  % RR3 IDs between 100 and 400
                    subjects = [subjects str2double(allDirs(jDir).name)];
                end
            end
            if isempty(subjects) 
                fprintf('NO VALID RR3 WM SUBJECTS FOUND\n')
                return
            end
        elseif iSubInput < 101 || iSubInput > 399  % RR3 IDs between 100 and 400
            clc
            fprintf('INVALID INPUT: %d\n', iSubInput)
            validSubId = true;  % stay in loop
            break
        end  % if iSubInput == 0
    end  % for iSubInput
end  % while validSubId

resolution = true; % Continues while loop until false
while resolution   % User input low vs high resolution
    plotResolution = 'Input [1] low resolution (10ms interval, 1Hz)\nInput [2] high resolution (2ms interval, 0.25Hz):';
    userPlot = input(plotResolution);
    resolution = false;  % exit loop
    if userPlot == 1
        timeArray    = lowResTimeArray;
        freqInterval = lowResFreqInterval;
        analysisType = '_low_res_';      % Becomes part of filename
        resLabel     = 'Low resolution'; % Displays when processing
        continue % Default is low resolution
    elseif userPlot == 2
        timeArray = highResTimeArray;
        freqInterval = highResFreqInterval;
        analysisType = '_high_res_';  % Becomes part of filename
        resLabel = 'High resolution'; % Displays when processing
    else % if invalid input
        userPlot = num2str(userPlot);
        fprintf('Invalid input: %s\n', userPlot)
        resolution = true;  % stay in loop
    end
end  % while resolution

tStart = tic; % Script start time
eeglab        % load eeglab

for iSub = 1:length(subjects) % loop through subjects
    tStart2     = tic; % Subject specific process time
    subId       = num2str(subjects(iSub)); % 3 digit sub id
    
    % Load input file
    fileName    = [subId, dataFileEnding];
    filePath    = [dataPath, subId, '/'];
    inputFile   = fullfile(filePath, fileName);
    
    if ~exist(inputFile, 'file')
        fprintf('MISSING %s INPUT FILE: %s\n', subId, inputFile)
        continue
    else
        EEG = pop_loadset('filename', fileName, 'filepath', filePath);
    end
    
    % Get item and order hit trials (exclude miss and no response trials)
    itemHitIndex  = [];
    orderHitIndex = [];
    itemMissIndex  = []; % by FH
    orderMissIndex = []; % by FH
    for jEpoch = 1:length(EEG.epoch) % loop through epochs
        for kLabel = 1:length(EEG.epoch(jEpoch).eventlabel) % loop through event labels (usually first label)
            labelName = EEG.epoch(jEpoch).eventlabel{kLabel};
            itemHitCheck  = regexp(labelName,'ITEM|hit');  % Creates 1x2 matrix if true
            orderHitCheck = regexp(labelName,'ORDER|hit'); % Creates 1x2 matrix if true
            itemMissCheck = regexp(labelName,'ITEM|miss');  % by FH
            orderMissCheck= regexp(labelName,'ORDER|miss'); % by FH
            if length(itemHitCheck) == 2
                itemHitIndex = [itemHitIndex jEpoch];
                break
            elseif length(orderHitCheck) == 2
                orderHitIndex  = [orderHitIndex jEpoch];
                break
            elseif length(itemMissCheck) == 2
                itemMissIndex= [itemMissIndex jEpoch]; % by FH
                break;
            elseif length(orderMissCheck) == 2
                orderMissIndex = [orderMissIndex jEpoch]; % by FH
                break;
            end
        end % kLabel
    end % jEpoch
    
    itemHitCount  = length(itemHitIndex);
    orderHitCount = length(orderHitIndex);
    itemMissCount  = length(itemMissIndex);   % by FH 
    orderMissCount = length(orderMissIndex);  % by FH
    
    % double check hit count (sanity check)
    if EEG.behaveItemHit ~= itemHitCount || EEG.behaveOrderHit ~= orderHitCount ...
            || EEG.behaveItemMiss ~= itemMissCount || EEG.behaveOrderMiss ~= orderMissCount % by FH
        fprintf('ERROR %s HIT COUNT: ITEM %d vs %d ORDER %d vs %d\n', subId, EEG.behaveItemHit, itemHitCount, EEG.behaveOrderHit, orderHitCount)
        fprintf('ERROR %s MISS COUNT: ITEM %d vs %d ORDER %d vs %d\n', subId, EEG.behaveItemMiss, itemMissCount, EEG.behaveOrderMiss, orderMissCount) % by FH
        continue
    else
        if itemHitCount > 0 && orderHitCount > 0 % See if hits exist
            %trialTypes={itemHitIndex orderHitIndex};
            trialTypes={itemHitIndex orderHitIndex itemMissIndex orderMissIndex};  % by FH
        else
            fprintf('MISSING HIT TRIALS: ITEM=%d ORDER=%d\n', itemHitCount, orderHitCount)
            continue
        end
    end
    
    % Output ERSP and ITC files
    erspFile = [filePath, subId, '_8_wm_ersp', analysisType, 'sr', num2str(tfCycles(1)), '_freqs_', num2str(freqRange(1)), '_', num2str(freqRange(2)), '_fh.mat'];
    itcFile  = [filePath, subId, '_8_wm_itc', analysisType, 'sr', num2str(tfCycles(1)), '_freqs_', num2str(freqRange(1)), '_', num2str(freqRange(2)), '_fh.mat'];
    
    % Header information for output files
    ERSP.subject       = subId;
    %ERSP.baseline      = tfBaseline;
    ERSP.baselineperiod = tfBaseline;
    %ERSP.conditions    = {'ITEM_HIT' 'ORDER_HIT' 'ORDER_v_ITEM_HIT'};
    ERSP.conditions    = {'ITEM_HIT' 'ORDER_HIT' 'ITEM_MISS' 'ORDER_MISS' 'ORDER_v_ITEM_HIT' 'ITEM_v_ORDER_HIT'}; % by FH. check 'trialTypes' variable
    ERSP.itemHitCount  = itemHitCount;
    ERSP.orderHitCount = orderHitCount;
    ERSP.itemMissCount  = itemMissCount;  % added by FH
    ERSP.orderMissCount = orderMissCount; % added by FH
    
    ITC.subject       = subId;
    %ITC.baseline      = tfBaseline;
    ITC.baselineperiod = tfBaseline;
    %ITC.conditions    = {'ITEM_HIT' 'ORDER_HIT' 'ORDER_v_ITEM_HIT'};
    ITC.conditions    = {'ITEM_HIT' 'ORDER_HIT' 'ITEM_MISS' 'ORDER_MISS' 'ORDER_v_ITEM_HIT' 'ITEM_v_ORDER_HIT'}; % by FH. check 'trialTypes' variable
    ITC.itemHitCount  = itemHitCount;
    ITC.orderHitCount = orderHitCount;
    ITC.itemMissCount  = itemMissCount;  % added by FH
    ITC.orderMissCount = orderMissCount; % added by FH
    
    for jElec = includedElectrodes % 64 electrodes (or manual edit)
        for kTrialType = 1:length(trialTypes) % Item hit or order hit
            fprintf('\nSUBJECT   : %s\nANALYSIS  : %s\nELECTRODE : %d/%d\nTRIAL TYPE: %s\n', subId, resLabel, jElec, length(includedElectrodes), ERSP.conditions{kTrialType})
            
% commented out by FH            
%             % eeglab newtimef function to create ERSPs
%             [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(EEG.data(jElec,:,trialTypes{kTrialType}), ...
%                 EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,tfCycles ,'freqs', freqRange, ...
%                 'elocs', EEG.chanlocs,'timesout',timeArray,'baseline',[NaN],'chaninfo', EEG.chaninfo,...
%                 'nfreqs',freqInterval,'plotphase','off','plotersp','off','plotitc','off');

            % added by FH
            if ~isempty(trialTypes{kTrialType})
                % eeglab newtimef function to create ERSPs
                [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(EEG.data(jElec,:,trialTypes{kTrialType}), ...
                    EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,tfCycles ,'freqs', freqRange, ...
                    'elocs', EEG.chanlocs,'timesout',timeArray,'baseline',[NaN],'chaninfo', EEG.chaninfo,...
                    'nfreqs',freqInterval,'plotphase','off','plotersp','off','plotitc','off');
            else
                tfdata = NaN;
            end
            % till here
            
            % Trial Data (baseline subtract on a trial by trial basis)
            [baselineVal1, baselineIndex1] = min(abs(times - tfBaseline(1))); % Value and index of first baseline point (or closest value to it)
            [baselineVal2, baselineIndex2] = min(abs(times - tfBaseline(2))); % Value and index of second baseline point (or closest value to it)
            baselineIndices                = baselineIndex1:baselineIndex2;   % Baseline time indices
            
            if baselineVal1 ~= 0 || baselineVal2 ~= 0 % Alert user if baseline values differ from input
                fprintf('\nCHANGED BASELINE FROM: %d:%d to %d:%d\n\n', tfBaseline(1), tfBaseline(2), times(baselineIndex1), times(baselineIndex2));
            end
            
% commented out by FH            
%             % Multiply complex data by conjugate and log transform to get ERSP for EACH trial (freqs, times, trials)
%             ERSP.noBaseline.cond{kTrialType}.chan{jElec} = double(10*log10(tfdata.*conj(tfdata)));
%             % Multiply complex data by conjugate and log transform to get ERSP for EACH trial baseline, then average baseline for all time points (freqs, 1, trials)
%             ERSP.subtractValues.cond{kTrialType}.chan{jElec} = mean(double(10*log10(tfdata(:,baselineIndices,:).*conj(tfdata(:,baselineIndices,:)))),2);
%             % Get data dimensions and pre-allocate 'ERSP.baseline' with [nan] values
%             dataDims = size(ERSP.noBaseline.cond{kTrialType}.chan{jElec});
%             ERSP.baseline.cond{kTrialType}.chan{jElec} = nan(dataDims(1), dataDims(2), dataDims(3));
%             
%             % Baseline subtract each trial with its own averaged baseline
%             for mTrial = 1:dataDims(3)
%                 for nCol = 1:dataDims(2)
%                     ERSP.baseline.cond{kTrialType}.chan{jElec}(:,nCol,mTrial) = ERSP.noBaseline.cond{kTrialType}.chan{jElec}(:,nCol,mTrial) - ERSP.subtractValues.cond{kTrialType}.chan{jElec}(:,1,mTrial);
%                 end
%             end
%             
%             % Average data across all trials
%             ERSP.noBaseline.cond{kTrialType}.chan{jElec}      = mean(ERSP.noBaseline.cond{kTrialType}.chan{jElec}, 3);
%             ERSP.subtractValues.cond{kTrialType}.chan{jElec}  = mean(ERSP.subtractValues.cond{kTrialType}.chan{jElec}, 3);
%             ERSP.baseline.cond{kTrialType}.chan{jElec}        = mean(ERSP.baseline.cond{kTrialType}.chan{jElec}, 3);
         
            % added by FH
            ERSP.noBaseline.cond{kTrialType}.chan{jElec}     = double(10*log10(mean(tfdata.*conj(tfdata),3)));
            ERSP.subtractValues.cond{kTrialType}.chan{jElec} = double(10*log10(mean(mean(tfdata(:,baselineIndices,:).*conj(tfdata(:,baselineIndices,:)),3),2)));
            for ifreq = 1:size(ERSP.noBaseline.cond{kTrialType}.chan{jElec},1)
                ERSP.baseline.cond{kTrialType}.chan{jElec}(ifreq,:)       = ERSP.noBaseline.cond{kTrialType}.chan{jElec}(ifreq,:) - ERSP.subtractValues.cond{kTrialType}.chan{jElec}(ifreq);
            end
            % till here
            
            ITC.trialTypes{kTrialType}.chan{jElec} = itc;  % Single trial itc    
            
        end  % kTrialType
% commented out by FH
%         % Order vs. Item Hit
%         ERSP.noBaseline.cond{3}.chan{jElec}     = ERSP.noBaseline.cond{2}.chan{jElec} - ERSP.noBaseline.cond{1}.chan{jElec};
%         ERSP.subtractValues.cond{3}.chan{jElec} = ERSP.subtractValues.cond{2}.chan{jElec} - ERSP.subtractValues.cond{1}.chan{jElec};
%         ERSP.baseline.cond{3}.chan{jElec}       = ERSP.baseline.cond{2}.chan{jElec} - ERSP.baseline.cond{1}.chan{jElec};
        
        % added by FH
        % Order vs. Item Hit and Item vs. Order Hit
        ERSP.noBaseline.cond{5}.chan{jElec}     = ERSP.noBaseline.cond{2}.chan{jElec} - ERSP.noBaseline.cond{1}.chan{jElec};
        ERSP.baseline.cond{5}.chan{jElec}       = ERSP.baseline.cond{2}.chan{jElec} - ERSP.baseline.cond{1}.chan{jElec};
        ERSP.noBaseline.cond{6}.chan{jElec}     = ERSP.noBaseline.cond{1}.chan{jElec} - ERSP.noBaseline.cond{2}.chan{jElec};
        ERSP.baseline.cond{6}.chan{jElec}       = ERSP.baseline.cond{1}.chan{jElec} - ERSP.baseline.cond{2}.chan{jElec};
        % till here
        
        ERSP.chanLabel{jElec} = EEG.chanlocs(jElec).labels; % electrode label
        ITC.chanLabel{jElec}  = EEG.chanlocs(jElec).labels; % electrode label
    end  % jElec
    
    ERSP.times       = times; % Time range used in analysis
    ERSP.freqs       = freqs; % Frequency range used in analysis
    ERSP.processTime = hourToc(tStart2, false); % Process time of subject
    
    ITC.times        = times; % Time range used in analysis
    ITC.freqs        = freqs; % Frequency range used in analysis
    ITC.processTime  = hourToc(tStart2, false); % Process time of subject
    
    % -v7.3 option saves files >2GB (not a problem with averaged data, but is a problem saving single trial data)
    save(erspFile, 'ERSP', '-v7.3');
    fprintf('CREATED: %s\n', erspFile);
    
    save(itcFile, 'ITC', '-v7.3');
    fprintf('CREATED: %s\nTIME ELAPSED: ', itcFile);
    
    hourToc(tStart2); % Display subject specific process time
    
end % iSub

fprintf('TOTAL TIME ELAPSED: ');
hourToc(tStart); % Display script process time
