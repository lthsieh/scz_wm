% RR3 WM ERSP GROUP (eeglab)
% Created: 06/24/2015 By: Evan Layher
% Revised: 11/13/2015 By: Evan Layher % Check time x freq dimensions
%
% Average subject time-frequency data per electrode per trial
% Outputs mat-file with control, patient and (control-patient) averages
% Input subject IDs below by 'controls' and 'patients'
%
% +++ NOTE: THIS SCRIPT USES A LOT OF MEMORY +++
% +++ WHEN TOTAL SUBJECT COUNT REACHES 80+ CONSIDERING MODIFYING +++

clear all
clc

% Control subjects to include (input numbers as 'strings')
controls = {'101' '102' '103' '104' ...
            '105' '106' '107' '108' ...
            '110' '111' '112' '113' ...
            '114' '115' '116' '117' ...
            '118' '119' '120' '121' ...
            '122' '123' '124' '125' ...
            '126' '127' '128'};

% Patient subjects to include (input numbers as 'strings')
patients = {'201' '202' '203' '204' ...
            '205' '206' '207' '208' ...
            '210' '211' '212' '213' ...
            '214' '215' '216' '218' ...
            '219' '220' '221' '222' ...
            '223'};

% Unmedicated patient subjects to include [put in 'patients' group if no separation of the 2 groups are desired]
noMeds = {'303' '304' '305' '306' ...
          '307' '308' '309' '310'};

freqInput = []; % [] = use frequencies in file, or include range (e.g. [3 40])
timeInput = []; % [] = use times in file, or include range (e.g. [-12800 3300])

wmDir    = '/nfs/to-eeg/wm_data/';
baseFile = '_8_wm_ersp_low_res_sr6_freqs_3_40_fh.mat'; % Files to load in
groupDir = [wmDir, 'group/'];
addpath('/nfs/to-eeg/code'); % Add code path

%--- 'inputData' index MUST correspond with 'outName' ---%
% inputData is the matlab structure name inside each 'baseFile' to average
inputDataTypes = {'ERSP.noBaseline' 'ERSP.baseline'};
inputFreqName = 'ERSP.freqs'; % Where frequencies are listed in 'inputDataTypes'
inputTimeName = 'ERSP.times'; % Where times are listed in 'inputDataTypes'

% outName are the output names of the averages
% NOTE: cn and sz counts and todays date are added to end of file
outNames = {'wm_no_baseline_sr6_freqs_3_40_' 'wm_baseline_sr6_freqs_3_40_'};

analysisDate = datestr(date, 'yyyy_mm_dd');
defaultIncludedElectrodes = 1:64; % Include all electrodes (1:64)

%%% ---Start of processing--- %%%
tStart = tic; % Record process time
if length(inputDataTypes) ~= length(outNames)
    error('Variables "inputData" and "outName" MUST BE SAME LENGTH\n')
end

allSubs = {controls{:} patients{:} noMeds{:}};

% Display amount of processing subjects
fprintf('COMBINING %d controls and %d patients...\n', length(controls), length(patients))
if ~isempty(noMeds) % Display if unmedicated patients specified
    fprintf('...and %d unmedicated patients...\n', length(noMeds))
end

if ~exist(groupDir, 'file')
    mkdir(groupDir)
end

for iSub = 1:length(allSubs) % Check that all data files exist
    subPath = [wmDir, allSubs{iSub}, '/', allSubs{iSub}, baseFile];
    if ~exist(subPath, 'file')
        error('MISSING ERSP FILE: %s\n', subPath) % Alert user of which files are missing
    end
end


lastCn = length(controls);           % End of control subjects
lastSz = lastCn + length(patients);  % End of patient subjects
lastNoMed = lastSz + length(noMeds); % End of unmedicated patients
for iInputType = 1:length(inputDataTypes) % baseline or no baseline
    subCount = 0; % Reset subject count
    inputData = inputDataTypes{iInputType};
    outName   = outNames{iInputType};
    
    if isempty(noMeds)
        outPath   = [groupDir, outName, num2str(length(controls)), 'CN_', num2str(length(patients)), 'SZ_', analysisDate, '_fh.mat'];
    else
        outPath   = [groupDir, outName, num2str(length(controls)), 'CN_', num2str(length(patients)), 'SZmed_', num2str(length(noMeds)), 'SZunmed_', analysisDate, '_fh.mat'];
    end
    
    fprintf('PROCESSING: %s ...\n', outPath)
    
    for jSub = 1:length(allSubs) % Loop through all subjects
        subCount = subCount + 1;
        subPath = [wmDir, allSubs{jSub}, '/', allSubs{jSub}, baseFile];
        fprintf('Loading data: %s\n', subPath) % Alert user that data is loading
        load(subPath); % Load subject data
        tempData = eval(inputData); % Baseline or No Baseline subject data
        tempFreqs = eval(inputFreqName); % Get frequency range from file
        tempTimes = eval(inputTimeName); % Get time range from file
        
        if isempty(tempData) || isempty(tempFreqs) || isempty(tempTimes)
            error('INVALID ARRAY NAME(s) IN "inputDataTypes", "inputFreqName" or "inputTimeName"\n')
        end
        
        % get frequency range
        if isempty(freqInput)
            finalFreq = tempFreqs; % All frequencies in file
            freqPlot  = 1:length(finalFreq); % frequency indices
        else % Use specified frequency range
            freqFind  = find(tempFreqs == freqInput(1) | tempFreqs == freqInput(2));
            if length(freqFind) == 2 % Must find 2 values, else input is not in file
                finalFreq = tempFreqs(freqFind(1):freqFind(2)); % Freq values
                freqPlot  = freqFind(1):freqFind(2); % Frequency indices
            else
                error('Invalid frequency range: %d - %d Hz Subject: %s\n', freqInput(1), freqInput(2), subPath)
            end
        end
        
        % get time range
        if isempty(timeInput)
            finalTime = tempTimes; % All file times
            timePlot = 1:length(finalTime); % Indices of time range
        else % Use specified time range
            timeFind = find(tempTimes == timeInput(1) | tempTimes == timeInput(2));
            if length(timeFind) == 2 % Must find 2 values, else input is not in file
                finalTime = tempTimes(timeFind(1):timeFind(2)); % time values
                timePlot  = timeFind(1):timeFind(2); % time indices
            else
                error('Invalid time range: %d - %d ms Subject: %s\n', timeInput(1), timeInput(2), subPath)
            end
        end
        
        if ~exist('combined', 'var') % Obtain values from first subject (should be same for all)
            if ~isempty(noMeds) % If unmedicated patients processed separately
                combined.patientsNoMeds = noMeds;
            end
            
            combined.chanLabels = ERSP.chanLabel;
            combined.conditions = ERSP.conditions;
            combined.controls   = controls; % List of controls in analysis
            combined.patients   = patients; % List of patients in analysis
            combined.freqs      = finalFreq;
            combined.times      = finalTime;
        end
        
        clear ERSP % Remove loaded ERSP file to save memory
        
        % Make sure frequencies and times match for each subject
        checkFreqs = intersect(combined.freqs, finalFreq);
        checkTimes = intersect(combined.times, finalTime);
        if length(checkFreqs) ~= length(finalFreq) || length(checkFreqs) ~= length(combined.freqs)
            error('Frequency mismatch %s\n', allSubs{jSub})            
        elseif length(checkTimes) ~= length(finalTime) || length(checkTimes) ~= length(combined.times)
            error('Time mismatch %s\n', allSubs{jSub})            
        end
        
        if jSub <= lastCn     % Control
            totalSubs = length(controls);
        elseif jSub <= lastSz % Patient
            totalSubs = length(patients);
        else                  % Unmedicated patient
            totalSubs = length(noMeds);
        end
        
        includedElectrodes = defaultIncludedElectrodes; % Reset electrode count for each subject
        includedElectrodes = removeElectrodes('wm', allSubs{jSub}, 1, includedElectrodes); % Do not include bad electrodes in average
        
        for kCon = 1:length(combined.conditions) % conditions: Item_hit, order_hit, item_miss, order_miss, order_v_item_hit item_v_order_hit
            for mChan = 1:length(combined.chanLabels) % electrodes: 64
                
                if subCount == 1 % Must create 'NaN' arrays for each electrode (to remove bad)
                    tempStruc{kCon}.chan{mChan} = nan(length(combined.freqs), length(combined.times), totalSubs); % Create NaN arrays to quicken processing
                end
                
                if ~isempty(intersect(includedElectrodes, mChan)) % Bad electrodes will be empty
                    tempStruc{kCon}.chan{mChan}(:,:,subCount) = tempData.cond{kCon}.chan{mChan}(freqPlot, timePlot);
                end
                
                if jSub == lastCn       % Control
                    if size(tempStruc{kCon}.chan{mChan}, 3) ~= legnth(controls)
                        error('Check Data Dimensions!');
                    end
                    combined.cnAvg.conditions{kCon}.chan{mChan} = nanmean(tempStruc{kCon}.chan{mChan}, 3); % Remove 'NaN' values (if bad electrodes exist)
                elseif jSub == lastSz   % Patient
                    if size(tempStruc{kCon}.chan{mChan}, 3) ~= legnth(patients)
                        error('Check Data Dimensions!');
                    end
                    combined.szAvg.conditions{kCon}.chan{mChan} = nanmean(tempStruc{kCon}.chan{mChan}, 3); % Remove 'NaN' values (if bad electrodes exist)
                    combined.cnVsz.conditions{kCon}.chan{mChan} = combined.cnAvg.conditions{kCon}.chan{mChan} - combined.szAvg.conditions{kCon}.chan{mChan};
                    tempStrucSZ{kCon}.chan{mChan} = tempStruc{kCon}.chan{mChan};
                elseif jSub == lastNoMed % Unmedicated patient
                    if size(tempStruc{kCon}.chan{mChan}, 3) ~= legnth(noMeds)
                        error('Check Data Dimensions!');
                    end
                    combined.sznomedsAvg.conditions{kCon}.chan{mChan} = nanmean(tempStruc{kCon}.chan{mChan}, 3); % Remove 'NaN' values (if bad electrodes exist)
                    combined.szsznomedsAvg.conditions{kCon}.chan{mChan} = nanmean(cat(3,tempStrucSZ{kCon}.chan{mChan},tempStruc{kCon}.chan{mChan}),3);
                    combined.cnVsznomeds.conditions{kCon}.chan{mChan} = combined.cnAvg.conditions{kCon}.chan{mChan} - combined.sznomedsAvg.conditions{kCon}.chan{mChan};
                    combined.cnVszszsznomeds.conditions{kCon}.chan{mChan} = combined.cnAvg.conditions{kCon}.chan{mChan} - combined.szsznomedsAvg.conditions{kCon}.chan{mChan};
                end
            end % mChan
        end % kCon
        
        if jSub == lastCn || jSub == lastSz || jSub == lastNoMed
            subCount = 0;   % Reset count
            clear tempStruc % clear temp structure
        end
    end % jSub
    
    fprintf('Saving...\n')
    save(outPath, 'combined', '-v7.3');
    fprintf('CREATED: %s\n', outPath)
    
    clear combined tempData % Remove structures to save memory
end % iInputType

hourToc(tStart);