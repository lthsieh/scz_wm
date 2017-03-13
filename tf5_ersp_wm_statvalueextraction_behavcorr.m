%% load data
clc
clear all
close all
path_file = '/Users/frank/Desktop/SCZ/';
load(fullfile(path_file,'wm_baseline_sr6_freqs_3_40_16CN_14SZmed_3SZunmed_2016_07_12_stats_fh.mat'));
%% find frequency and time ranges
start_time = [0 3000]; % in ms, relative to onset of delay fixation cross
start_freq = [8 10];    % frequency range
index_time = [find(precombined.times == start_time(1)) find(precombined.times == start_time(2))];
index_freq = [find(precombined.freqs == start_freq(1)) find(precombined.freqs == start_freq(2))];
if length(index_time) ~= 2 || length(index_freq) ~= 2
    error('some frequencies or time points could not be found!');
end

%% conditions of interest
conditions_to_extract = {'ITEM_HIT','ITEM_MISS','ORDER_HIT','ORDER_MISS'};
condition_index = [];
for icon = 1:length(conditions_to_extract)
    tmp = find(strcmp(conditions_to_extract{icon},precombined.conditions));
    condition_index = [condition_index tmp];
end
if length(condition_index) ~= length(conditions_to_extract)
    error('Some conditions were not found!');
end

%% participant groups
participant_groups = {'cn','sz','sznomeds','szsznomeds'};

%% Group electrodes (Variable names will become titles in graphs)
groupRegions = {'Left_PFC',     [1 4 6 7 8 15 16 17];...
                'Mid_PFC',      [2 9 10 11 18 19 20];...
                'Right_PFC',    [3 5 12 13 14 21 22 23];...
                'Left_Central', [24 25 26 34 35 36];...
                'Mid_Central',  [27 28 29 37 38 39];...
                'Right_Central',[30 31 32 40 41 42];...
                'Left_Post',    [44 45 46 53 54];...
                'Mid_Post',     [47 48 49 55 56 57 61 62 63];...
                'Right_Post',   [50 51 52 58 59]};

%%
for igroup = 1:length(participant_groups)
    switch participant_groups{igroup}
        case 'cn'
            data = precombined.cn;
        case 'sz'
            data = precombined.sz;
        case 'sznomeds'
            data = precombined.sznomeds;
        case 'szsznomeds'
            data = precombined.szsznomeds;
    end
    clear data_excel;
    column_counter = 0;
    for icon = condition_index
        column_counter = column_counter + 1;
        for iregion = 1:size(groupRegions,1)
            clear data_region;
            index_channels = groupRegions{iregion,2};
            for ichan = 1:length(index_channels)
                data_region(:,:,:,ichan) = data.conditions{icon}.chan{index_channels(ichan)};
            end
            data_region = nanmean(data_region,4);
            data_region = data_region(index_freq(1):index_freq(2),index_time(1):index_time(2),:); % extract desired freq and time range
            data_region = nanmean(data_region,1); % average across freqs
            data_region = nanmean(data_region,2); % average across times
            data_region = squeeze(data_region); % get averaged time-frequency data for individual participants
            data_excel(:,(iregion-1)*length(condition_index) + column_counter) = data_region;
        end
    end
    switch participant_groups{igroup}
        case 'cn'
            data_excel_cn = data_excel;
        case 'sz'
            data_excel_sz = data_excel;
        case 'sznomeds'
            data_excel_sznomeds = data_excel;
        case 'szsznomeds'
            data_excel_szsznomeds = data_excel;
    end
end

%% power-behavioral correlations


