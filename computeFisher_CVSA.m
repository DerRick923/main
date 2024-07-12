% compute fischer score matrix with dimensions [freq_bands x channels x runs]
close all
clear all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA\cnbi-smrtrain\toolboxes\cva'))

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
       '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
 
% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'POZ', 'O1', 'O2', '', ...
%        '', '', '', '', '', '', '', '', '', '', '', '', '', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
% file info
subject = 'c7';
lap_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Laplacian\lap_39ch_CVSA.mat';
chanlocs_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Chanlocs\new_chanlocs64.mat';
path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' subject '\mat_selectedTrials'];
% path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\ ' subject '\gdf'];
% ubuntu path = ['/home/riccardo/test_ws/records/ ' subject '/matselected_Trials];
classLb = {'Task 1','Task 2'};
classes = [730,731];

matfiles = dir(fullfile(path, '*.mat'));
load(lap_path);
load(chanlocs_path);

band = {[6 9], [9 12], [8 14], [12 15], [15 18], [18 21]};
nbands = length(band);

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]); Rk=[];
for i=1:length(matfiles)
    file = fullfile(path, matfiles(i).name);
    load(file);
    curr_s = signal(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali
    slap = curr_s*lap;
    curr_s = slap;
    curr_h = header.EVENT;
    %Create Rk vector (run)
    cRk = i*ones(size(curr_s,1),1);
    Rk = cat(1,Rk,cRk);
    % concateno eventi
    events.TYP = cat(1, events.TYP, curr_h.TYP);
    events.DUR = cat(1, events.DUR, curr_h.DUR);
    events.POS = cat(1, events.POS, curr_h.POS + size(s, 1));
    s = cat(1, s, curr_s);
end

%% Create Vector labels
[nsamples,nchannels] = size(s);
[feedb_pos, feedb_dur, fix_dur, fix_pos, cue_dur, cue_pos, ntrials] = extract_info_label(events, 781, 786, [730 731]);

%% Extract trial data
[TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(s, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, ntrials);

%% Data processing
s_processed = NaN(nsamples,nchannels,nbands);
for f_idx=1:nbands
    sel_band = band{f_idx}; %Hz
    t_window = 1; %[s]
    windowSize = events.SampleRate*t_window;
    filtOrder = 4;
    s_movavg = data_processing(s, nchannels, events.SampleRate, sel_band, filtOrder, t_window);
    s_processed(:,:,f_idx) = s_movavg;
end
%% Trial extraction
%si estraggono i dati dalla fixation alla fine del feedback
trial_dur = min(TrialStop-TrialStart);
TrialData = []; new_Rk = []; new_Ck = [];
tCk = zeros(ntrials,1);
for trId=1:ntrials
    cstart = TrialStart(trId);
    cstop = cstart + trial_dur - 1;

%TrialData(:,:,:,trId) = s_processed(cstart:cstop,:,:);

    c_TrialData = s_processed(cstart:cstop,:,:);
    TrialData = cat(1,TrialData,c_TrialData);
    c_Rk = Rk(cstart:cstop,1);
    new_Rk = cat(1,new_Rk,c_Rk);
    c_Ck = Ck(cstart:cstop);
    new_Ck = cat(1,new_Ck,c_Ck);
    tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
end

%% Baseline extraction
% minFix_dur = min(FixStop - FixStart);
% Reference = NaN(minFix_dur, nchannels, ntrials);
% for trId=1:ntrials
%     cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
%     cstop = cstart+ minFix_dur - 1;
%     Reference(:,:,trId) = movavg_alpha(cstart:cstop,:);
% end
%Baseline = repmat(mean(Reference),[size(TrialData,1) 1 1]);

%% Compute ERD and LogBandPOwer [samples x channels x bands] con tutti i trial
%ERD = log(TrialData./Baseline);
ERD = log(TrialData);       %Logband
ERD = permute(ERD, [1 3 2]);

%% Compute Fischer score
Runs = unique(Rk);
nruns = length(Runs);

classes = [730, 731];
nclasses = length(classes);
% fischer_score = NaN(nbands,nchannels,nruns);
% F2S = NaN(nbands*nchannels,nruns);
% for rId=1:nruns
%     rindex = new_Rk==Runs(rId);
%     cmu = NaN(nbands,nchannels,2);
%     csigma = NaN(nbands,nchannels,2);
% 
%     for cId=1:nclasses
%        cindex = rindex & new_Ck==classes(cId);
%        cmu(:,:,cId) = squeeze(mean(ERD(cindex,:,:)));
%        csigma(:,:,cId) = squeeze(std(ERD(cindex,:,:)));
%     end
%     fischer_score(:,:,rId) = abs(cmu(:,:,2)-cmu(:,:,1))./sqrt((csigma(:,:,1).^2 + csigma(:,:,2).^2));
% end
%% Compute CVA
cva = nan(nbands, nchannels, nruns);
%cva = nan(length(intervals)+2, length(idx_selFreqs), nchannels);
for idx_r = 1:nruns
    for i= 1:nbands
        rindex = new_Rk==Runs(idx_r);
        c_data = squeeze(ERD(rindex,i,:));
        c_ck = new_Ck(rindex);
        c = cva_tun_opt(c_data, c_ck);
        cva(i, :, idx_r) = c;
    end
end

%% Visualization
%% Visualization Fisher score
disp('[proc] |- Visualizing fisher score for offline runs');
freq_intervals = {'6-9', '9-12', '8-14', '12-15', '15-18', '18-21'};
OfflineRuns = unique(new_Rk);
NumCols = length(OfflineRuns);
climits = [];
handles = nan(length(OfflineRuns), 1);
a = find(~strcmp(channels_label,''));
fig1 = figure;
colormap('jet');
for rId = 1:length(OfflineRuns)
        subplot(2, ceil(NumCols/2), rId);
        imagesc(cva(:, a, OfflineRuns(rId))');
        axis square;
        colorbar;
        set(gca, 'XTick', 1:nbands);
        set(gca, 'XTickLabel', freq_intervals);
        set(gca, 'YTick', 1:size(a,2));
        set(gca, 'YTickLabel', channels_label(find(~strcmp(channels_label,''))));
        xtickangle(90);
        xlabel('Hz');
        ylabel('channel');
        
        title(['Calibration run ' num2str(OfflineRuns(rId))]);
        
        climits = cat(2, climits, get(gca, 'CLim'));
        handles(OfflineRuns(rId)) = gca;
    end
    
    
        set(handles, 'clim', [0 max(max(climits))]);
        sgtitle(['Fisher score Subj: ' subject]);

% fisher_score_total = NaN(nbands,nchannels);
    % %F2S_total = NaN(nbands*nchannels,nruns);
    % cmu_total = NaN(nbands,nchannels,2);
    % csigma_total = NaN(nbands,nchannels,2);
    % for cId=1:nclasses
    %        cindex_new = new_Ck==classes(cId);
    %        cmu_total(:,:,cId) = squeeze(mean(ERD(cindex_new,:,:)));
    %        csigma_total(:,:,cId) = squeeze(std(ERD(cindex_new,:,:)));
    % end
    % fischer_score_total(:,:) = abs(cmu_total(:,:,2)-cmu_total(:,:,1))./sqrt((csigma_total(:,:,1).^2 + csigma_total(:,:,2).^2));
cva_total = NaN(nbands,nchannels);
for i= 1:nbands
        c_data = squeeze(ERD(:,i,:));
        c_ck = new_Ck;
        c = cva_tun_opt(c_data, c_ck);
        cva_total(i, :) = c;
end
fig2=figure;
colormap('jet');
imagesc(cva_total(:,a)');
axis square;
colorbar;
set(gca, 'XTick', 1:nbands);
set(gca, 'XTickLabel', freq_intervals);
set(gca, 'YTick', 1:size(a,2));
set(gca, 'YTickLabel', channels_label(find(~strcmp(channels_label,''))));
xtickangle(90);
xlabel('Hz');
ylabel('channel');
title(['Total FS Subj: ' subject]);



% saving fischer score for UI
rowLabels = channels_label(find(~strcmp(channels_label,'')));
colLabels = freq_intervals;
cva_selected = cva_total(:,a)';
save('fischer_scores.mat', 'cva_selected', 'rowLabels', 'colLabels');
   
