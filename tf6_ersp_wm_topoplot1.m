%% load data
clc
clear all
close all
path_file = '/Users/frank/Desktop/SCZ/';
BASELINE_CORRECTION = true;
if BASELINE_CORRECTION
    load(fullfile(path_file,'wm_baseline_sr6_freqs_3_40_27CN_21SZmed_8SZunmed_2017_02_15_stats_fh.mat'));
else
    load(fullfile(path_file,'wm_no_baseline_sr6_freqs_3_40_27CN_21SZmed_8SZunmed_2017_02_15_stats_fh.mat'));
end

%% participant groups
participant_groups = {'cn','sz','sznomeds','szsznomeds'};

%% topomaps
TMAPS_OR_POWERVALUE = 'tmaps'; % 'powervalue'

%% find frequency and time ranges
start_time = [0 3000]; % in ms, relative to onset of delay fixation cross
start_freq = [4 8];    % frequency range
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

%% obtain values for individual conditions from individual participants
for igroup = 1:length(participant_groups)
    switch participant_groups{igroup}
        case 'cn'
            data = precombined.cn;       % precombined.cn.conditions{1}.chan{1}: frequency x time x participant
        case 'sz'
            data = precombined.sz;       % precombined.sz.conditions{1}.chan{1}: frequency x time x participant
        case 'sznomeds'
            data = precombined.sznomeds; % precombined.sznomeds.conditions{1}.chan{1}: frequency x time x participant
        case 'szsznomeds'
            data = precombined.szsznomeds; % precombined.sznomeds.conditions{1}.chan{1}: frequency x time x participant
    end
    
    clear data_topo;
    for icon = 1:length(conditions_to_extract)
        data_topo{1,icon} = conditions_to_extract{icon};
        chan_num = length(data.conditions{condition_index(icon)}.chan);
        for ichan = 1:chan_num
            tmp = data.conditions{condition_index(icon)}.chan{ichan}; % freq x time x participants
            tmp = tmp(index_freq(1):index_freq(2),index_time(1):index_time(2),:); % extract desired freq and time range
            tmp = nanmean(tmp,1); % average across freqs
            tmp = nanmean(tmp,2); % average across times
            tmp = squeeze(tmp);   % get individual participants' values
            data_topo{2,icon}(:,ichan) = tmp;
        end
    end
    switch participant_groups{igroup}
        case 'cn'
            data_topo_cn       = data_topo;
        case 'sz'
            data_topo_sz       = data_topo;
        case 'sznomeds'
            data_topo_sznomeds = data_topo;
        case 'szsznomeds'
            data_topo_szsznomeds = data_topo;
    end
end

%% contrasts - within individual groups
contrasts = {{'ITEM_HIT' ,'ORDER_HIT'};...  % the left minus the right
             {'ORDER_HIT','ITEM_HIT'}};
         
clear topo_cn topo_sz topo_sznomeds topo_szsznomeds
for igroup = 1:4
    clear topo_tmp;
    if igroup == 1
        data_topo = data_topo_cn;
    elseif igroup == 2
        data_topo = data_topo_sz;
    elseif igroup == 3
        data_topo = data_topo_sznomeds;
    elseif igroup == 4
        data_topo = data_topo_szsznomeds;
    end
    switch TMAPS_OR_POWERVALUE
     case 'tmaps'
         for icontrast = 1:length(contrasts)
             topo_tmp{1,icontrast} = [contrasts{icontrast,1}{1},'-',contrasts{icontrast,1}{2}];
             index1 = find(strcmp(contrasts{icontrast}{1},data_topo(1,:)));
             index2 = find(strcmp(contrasts{icontrast}{2},data_topo(1,:)));
             tmp_data1 = data_topo{2,index1};
             tmp_data2 = data_topo{2,index2};
             for ichan = 1:size(tmp_data1,2)
                 [h,p,ci,stat] = ttest(tmp_data1(:,ichan),tmp_data2(:,ichan));
                 topo_tmp{2,icontrast}(1,ichan) = stat.tstat;
             end
         end
     case 'powervalue'
         for icontrast = 1:length(contrasts)
             topo_tmp{1,icontrast} = [contrasts{icontrast,1}{1},'-',contrasts{icontrast,1}{2}];
             index1 = find(strcmp(contrasts{icontrast}{1},data_topo(1,:)));
             index2 = find(strcmp(contrasts{icontrast}{2},data_topo(1,:)));
             tmp_data1 = nanmean(data_topo{2,index1},1);
             tmp_data2 = nanmean(data_topo{2,index2},1);
             for ichan = 1:size(tmp_data1,2)
                 topo_tmp{2,icontrast}(1,ichan) = tmp_data1(1,ichan) - tmp_data2(1,ichan);
             end
         end
    end
    if igroup == 1
        topo_cn          = topo_tmp;
    elseif igroup == 2
        topo_sz          = topo_tmp;
    elseif igroup == 3
        topo_sznomeds    = topo_tmp;
    elseif igroup == 4
        topo_szsznomeds = topo_tmp;
    end
end                  

%% topoplot
group = 'SZ'; %'CN' 'SZ' 'SZNOMEDS' 'SZNSZNOMEDS'
switch group
    case 'CN'
        tmp_data = topo_cn;
    case 'SZ'
        tmp_data = topo_sz;
    case 'SZNOMEDS'
        tmp_data = topo_sznomeds;
    case 'SZNSZNOMEDS'
        tmp_data = topo_szsznomeds;
end
datavector = tmp_data{2,2};
datatitle  = tmp_data{1,2};
figure;
topoplot(datavector, EEG.chanlocs);

title([group,' ',strrep(datatitle,'_',' '),' ',[num2str(start_freq(1)),'-',num2str(start_freq(2)),'Hz'],' ',[num2str(start_time(1)),'-',num2str(start_time(2)),'ms']]);
colorbar;


%% contrasts - between groups



