%% load in data
%% testtest
clc
clear all
close all
filepath = '/Volumes/LaCie/SCZ/';
filename = 'wm_baseline_sr6_freqs_3_40_27CN_21SZmed_8SZunmed_2017_03_01_stats_fh.mat';
load(fullfile(filepath,filename));

%% organize data for hierarchical clustering
% controls
data_control = [];
for ichan = 1:length(precombined.chanLabels)
    clear tmp_data;
    for icon = 1:length(precombined.conditions)
        tmp_data(:,:,icon) = nanmean(precombined.cn.conditions{icon}.chan{ichan},3); % average across controls
    end
    tmp_data = mean(tmp_data,3);                % average across conditions
    data_control = [data_control,tmp_data];     % concatenate time-frequency data across electrodes
end
% sz & sznomeds
data_szsznomeds = [];
for ichan = 1:length(precombined.chanLabels)
    clear tmp_data;
    for icon = 1:length(precombined.conditions)
        tmp = cat(3,precombined.sz.conditions{icon}.chan{ichan},precombined.sznomeds.conditions{icon}.chan{ichan});
        tmp_data(:,:,icon) = nanmean(tmp,3);       % average across patients
    end
    tmp_data = mean(tmp_data,3);                % average across conditions
    data_szsznomeds = [data_szsznomeds,tmp_data];     % concatenate time-frequency data across electrodes
end
% sz
data_sz = [];
for ichan = 1:length(precombined.chanLabels)
    clear tmp_data;
    for icon = 1:length(precombined.conditions)
        tmp_data(:,:,icon) = nanmean(precombined.sz.conditions{icon}.chan{ichan},3); % average across patients
    end
    tmp_data = mean(tmp_data,3);                % average across conditions
    data_sz = [data_sz,tmp_data];     % concatenate time-frequency data across electrodes
end
% sznomeds
data_sznomeds = [];
for ichan = 1:length(precombined.chanLabels)
    clear tmp_data;
    for icon = 1:length(precombined.conditions)
        tmp_data(:,:,icon) = nanmean(precombined.sznomeds.conditions{icon}.chan{ichan},3); % average across patients
    end
    tmp_data = mean(tmp_data,3);                % average across conditions
    data_sznomeds = [data_sznomeds,tmp_data];     % concatenate time-frequency data across electrodes
end

%% 
frequencies = precombined.freqs';
fontsize = 16;

%% hierarchical clustering over frequencies - controls
data_control_cluster = linkage(data_control,'average','correlation');  % equivalent to ### euc_dist = pdist(data_control,'correlation'); clust_tree_euc = linkage(euc_dist, 'average'); ###
data_control_pdist = pdist(data_control,'correlation'); % default is Euclidean distance
leafOrder = optimalleaforder(data_control_cluster,data_control_pdist);
% plot dendrogram
figure; 
[control_H,control_T,control_outperm]=dendrogram(data_control_cluster,'Reorder',leafOrder,'ColorThreshold','default');
title('Controls','Fontsize',fontsize)
set(control_H,'LineWidth',2);
h1 = gca;
set(h1,'Ticklength',[0 0]);                  % remove tick marks
old_XTicklabel = str2num(h1.XTickLabel);     % get XTick labels
clear new_XTicklabel
for i = 1:length(old_XTicklabel)
    new_XTicklabel{i} = mean(frequencies(find(control_T == old_XTicklabel(i))));
end
set(h1,'XTickLabel',new_XTicklabel,'XTickLabelRotation',90);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Distance','Fontsize',fontsize);
% plot dissimilarity matrix
data_control_pdist_square = squareform(data_control_pdist);
figure;
imagesc(data_control_pdist_square);
h1 = gca;
old_XTicklabel = str2num(str2mat(h1.XTickLabel{:}));     % get XTick labels
new_XTicklabel = old_XTicklabel + 2;
new_YTicklabel = new_XTicklabel;
set(gca,'XTickLabel',num2cell(new_XTicklabel),'YTickLabel',num2cell(new_YTicklabel));
title('Controls - Dissimilarity Matrix','Fontsize',fontsize);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Frequency (Hz)','Fontsize',fontsize);
caxis([0 1]);
colorbar;


%% hierarchical clustering over frequencies - szsznomeds
data_szsznomeds_cluster = linkage(data_szsznomeds,'average','correlation');
data_szsznomeds_pdist = pdist(data_szsznomeds,'correlation'); % default is Euclidean distance
leafOrder = optimalleaforder(data_szsznomeds_cluster,data_szsznomeds_pdist);
% plot dendrogram
figure; 
[patient_H,patient_T,patient_outperm]=dendrogram(data_szsznomeds_cluster,'Reorder',leafOrder,'ColorThreshold','default'); 
title('Patients','Fontsize',fontsize)
set(patient_H,'LineWidth',2);
h2 = gca;
set(h2,'Ticklength',[0 0]);                  % remove tick marks
old_XTicklabel = str2num(h2.XTickLabel);     % get XTick labels
clear new_XTicklabel
for i = 1:length(old_XTicklabel)
    new_XTicklabel{i} = mean(frequencies(find(patient_T == old_XTicklabel(i))));
end
set(h2,'XTickLabel',new_XTicklabel,'XTickLabelRotation',90);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Distance','Fontsize',fontsize);
% plot dissimilarity matrix
data_szsznomeds_pdist_square = squareform(data_szsznomeds_pdist);
figure;
imagesc(data_szsznomeds_pdist_square);
h2 = gca;
old_XTicklabel = str2num(str2mat(h2.XTickLabel{:}));     % get XTick labels
new_XTicklabel = old_XTicklabel + 2;
new_YTicklabel = new_XTicklabel;
set(gca,'XTickLabel',num2cell(new_XTicklabel),'YTickLabel',num2cell(new_YTicklabel));
title('Patients - Dissimilarity Matrix','Fontsize',fontsize);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Frequency (Hz)','Fontsize',fontsize);
caxis([0 1]);
colorbar;


%% hierarchical clustering over frequencies - sz
data_sz_cluster = linkage(data_sz,'average','correlation');
data_sz_pdist = pdist(data_sz,'correlation'); % default is Euclidean distance
leafOrder = optimalleaforder(data_sz_cluster,data_sz_pdist);
% plot dendrogram
figure; 
[patient_H,patient_T,patient_outperm]=dendrogram(data_sz_cluster,'Reorder',leafOrder,'ColorThreshold','default'); 
title('Patients','Fontsize',fontsize)
set(patient_H,'LineWidth',2);
h2 = gca;
set(h2,'Ticklength',[0 0]);                  % remove tick marks
old_XTicklabel = str2num(h2.XTickLabel);     % get XTick labels
clear new_XTicklabel
for i = 1:length(old_XTicklabel)
    new_XTicklabel{i} = mean(frequencies(find(patient_T == old_XTicklabel(i))));
end
set(h2,'XTickLabel',new_XTicklabel,'XTickLabelRotation',90);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Distance','Fontsize',fontsize);
% plot dissimilarity matrix
data_sz_pdist_square = squareform(data_sz_pdist);
figure;
imagesc(data_sz_pdist_square);
h2 = gca;
old_XTicklabel = str2num(str2mat(h2.XTickLabel{:}));     % get XTick labels
new_XTicklabel = old_XTicklabel + 2;
new_YTicklabel = new_XTicklabel;
set(gca,'XTickLabel',num2cell(new_XTicklabel),'YTickLabel',num2cell(new_YTicklabel));
title('Patients - Dissimilarity Matrix','Fontsize',fontsize);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Frequency (Hz)','Fontsize',fontsize);
caxis([0 1]);
colorbar;


%% hierarchical clustering over frequencies - sznomeds
data_sznomeds_cluster = linkage(data_sznomeds,'average','correlation');
data_sznomeds_pdist = pdist(data_sznomeds,'correlation'); % default is Euclidean distance
leafOrder = optimalleaforder(data_sznomeds_cluster,data_sznomeds_pdist);
% plot dendrogram
figure; 
[patient_H,patient_T,patient_outperm]=dendrogram(data_sznomeds_cluster,'Reorder',leafOrder,'ColorThreshold','default'); 
title('Patients','Fontsize',fontsize)
set(patient_H,'LineWidth',2);
h2 = gca;
set(h2,'Ticklength',[0 0]);                  % remove tick marks
old_XTicklabel = str2num(h2.XTickLabel);     % get XTick labels
clear new_XTicklabel
for i = 1:length(old_XTicklabel)
    new_XTicklabel{i} = mean(frequencies(find(patient_T == old_XTicklabel(i))));
end
set(h2,'XTickLabel',new_XTicklabel,'XTickLabelRotation',90);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Distance','Fontsize',fontsize);
% plot dissimilarity matrix
data_sznomeds_pdist_square = squareform(data_sznomeds_pdist);
figure;
imagesc(data_sznomeds_pdist_square);
h2 = gca;
old_XTicklabel = str2num(str2mat(h2.XTickLabel{:}));     % get XTick labels
new_XTicklabel = old_XTicklabel + 2;
new_YTicklabel = new_XTicklabel;
set(gca,'XTickLabel',num2cell(new_XTicklabel),'YTickLabel',num2cell(new_YTicklabel));
title('Patients - Dissimilarity Matrix','Fontsize',fontsize);
xlabel('Frequency (Hz)','Fontsize',fontsize);
ylabel('Frequency (Hz)','Fontsize',fontsize);
caxis([0 1]);
colorbar;


%% dendrogram - controls
euc_dist = pdist(data_control,'correlation');
clust_tree_euc = linkage(euc_dist, 'average');
figure;
dendrogram(clust_tree_euc);


%% dendrogram - szsznomeds
euc_dist = pdist(data_szsznomeds,'correlation');
clust_tree_euc = linkage(euc_dist, 'average');
figure;
dendrogram(clust_tree_euc);

%% dendrogram - sz
euc_dist = pdist(data_sz,'correlation');
clust_tree_euc = linkage(euc_dist, 'average');
figure;
dendrogram(clust_tree_euc);

%% dendrogram - sznomeds
euc_dist = pdist(data_sznomeds,'correlation');
clust_tree_euc = linkage(euc_dist, 'average');
figure;
dendrogram(clust_tree_euc);






