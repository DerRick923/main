% compute fischer score matrix with dimensions [freq_bands x channels x runs]
close all
clear all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\eeglab2024.0'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))
%addpath(genpath('/media/riccardo/A658ED4B58ED1B37/Users/User/Desktop/MATLAB/CVSA'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA\cnbi-smrtrain\toolboxes\cva'))
%addpath(genpath('/media/riccardo/A658ED4B58ED1B37/Users/User/Desktop/MATLAB/CVSA/cnbi-smrtrain/toolboxes/cva'))

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
       '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
 
% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};


% file info
c_subject = 'h7';
lap_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Laplacian\lap_39ch_CVSA.mat';
chanlocs_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Chanlocs\new_chanlocs64.mat';
path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' c_subject '\mat_selectedTrials'];
%path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' c_subject '\gdf'];

% ubuntu path = ['/home/riccardo/test_ws/records/ ' c_subject '/matselected_Trials];
%lap_path = '/media/riccardo/A658ED4B58ED1B37/Users/User/Desktop/MATLAB/CVSA/Laplacian/lap_39ch_CVSA.mat';
%chanlocs_path = '/media/riccardo/A658ED4B58ED1B37/Users/User/Desktop/MATLAB/CVSA/Chanlocs/new_chanlocs64.mat';

classLb = {'Task 1','Task 2'};
classes = [730,731];
nclasses = length(classes);

load(lap_path);
load(chanlocs_path);

%files = dir(fullfile(path, '*.gdf'));  %for ubuntu and gdf
files = dir(fullfile(path, '*.mat'));

band = {[6 9], [9 12], [8 14], [12 15], [15 18], [18 21]};
nbands = length(band);

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]); Rk=[];
for i=1:length(files)
    file = fullfile(path, files(i).name);

    load(file);
    %[signal,header] = sload(file);
    
    curr_s = signal(:,1:39);
    %slap = curr_s*lap;
    %curr_s = slap;
    curr_h = header.EVENT;
    isgdf = contains(file, '.gdf');
    iscalibration = contains(file, 'calibration');
    if isgdf && iscalibration
        start = find(curr_h.TYP == 1,1,'first');
        curr_h.TYP = curr_h.TYP(start:end);
        curr_h.POS = curr_h.POS(start:end);
        curr_h.DUR = curr_h.DUR(start:end);
    end
    % Create Rk vector (run)
    cRk = i*ones(size(curr_s,1),1);
    Rk = cat(1,Rk,cRk);
    % Concatenate events
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
DataperTrial = NaN(trial_dur,nchannels,nbands,ntrials);
tCk = zeros(ntrials,1);
for trId=1:ntrials
    cstart = TrialStart(trId);
    cstop = cstart + trial_dur - 1;

    DataperTrial(:,:,:,trId) = s_processed(cstart:cstop,:,:);

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
ERD_p= permute(ERD, [1 3 2]);

%% Compute Fischer score
Runs = unique(Rk);
nruns = length(Runs);

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
        c_data = squeeze(ERD_p(rindex,i,:));
        c_ck = new_Ck(rindex);
        c = cva_tun_opt(c_data, c_ck);
        cva(i, :, idx_r) = c;
    end
end

%% Visualization
%% Visualization Fisher score
disp('[proc] |- Visualizing fischer score for offline runs');
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
sgtitle(['Fisher score Subj: ' c_subject]);

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
        c_data = squeeze(ERD_p(:,i,:));
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
title(['Total FS Subj: ' c_subject]);



%% Topoplot LogBand power
ERDpertrial = log(DataperTrial);
chanlocs_label = {chanlocs.labels};
fixPeriod = [1/events.SampleRate 2]*events.SampleRate;
cuePeriod = [2 3]*events.SampleRate;
cfPeriod = [3 4; 4 5; 5 6; 6 (trial_dur/events.SampleRate)]*events.SampleRate;
    
%Select channels
recorded_channels =  {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
'', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
    cuetocf_erd1 = mean(mean(ERDpertrial(cuePeriod(1):cfPeriod(2), :, :, tCk == classes(1)), 4), 1);
    cuetocf_erd2 = mean(mean(ERDpertrial(cuePeriod(1):cfPeriod(2), :, :, tCk == classes(2)), 4), 1);
    cuetocf_erd = squeeze(cuetocf_erd2-cuetocf_erd1);
    cuetoc_feed = zeros(64, nbands);
        for i=1:length(chanlocs_label)
            for j = 1:length(recorded_channels)
                if strcmpi(chanlocs_label{i}, recorded_channels{j})
                    if ~isnan(cuetocf_erd(j))
                    cuetoc_feed(i,:) = cuetocf_erd(j,:);
                    else
                    cuetoc_feed(i,:) = 0;
                    end
                end
            end
        end
% cuetocfeed contiene i valori della logband power per ogni frequency band
% corrispondenti ai canali che vengono registrati

% saving subject id
subj = 'C:\Users\User\Desktop\MATLAB\CVSA\main\c_subject.mat';
save(subj,'c_subject');

% saving fischer score for UI
feature_file = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' c_subject '\dataset\fischer_scores.mat'];
% feature_file = [/home/riccardo/test_ws/records/' c_subject '/dataset/fischer_scores.mat'];
rowLabels = channels_label(find(~strcmp(channels_label,'')));
colLabels = freq_intervals;
cva_selected = cva_total(:,a)';
save(feature_file, 'cva_selected', 'rowLabels', 'colLabels','band');

% saving logband for UI
logband_file = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' c_subject '\dataset\logband_power.mat'];
% logband_file = [/home/riccardo/test_ws/records/' c-subject '/dataset/logband_power.mat'];
logbandPower = cuetoc_feed;
electrodePos = chanlocs;
save(logband_file,'logbandPower','electrodePos');