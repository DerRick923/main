close all
clear all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))

%channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
%        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
        'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% file info
subject = 'g2';
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

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]);
band = {[6 9], [9 12], [8 14], [12 15], [15 18], [18 21]};

for f_idx=1:length(band)
    sel_band = band{f_idx}; %[hz]
    for i=1:length(matfiles)
        file = fullfile(path, matfiles(i).name);
        load(file);
        curr_s = signal(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali
        slap = curr_s*lap;
        curr_s = slap;
        curr_h = header.EVENT;
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
    t_window = 1; %[s]
    windowSize = events.SampleRate*t_window;
    % band = [8 14]; %[hz]
    filtOrder = 4;
    movavg_alpha = data_processing(s,nchannels,events.SampleRate,sel_band,filtOrder,t_window);
    
    %% Trial extraction
    trial_dur = min(TrialStop-TrialStart);
    TrialData = NaN(trial_dur,nchannels,ntrials);
    tCk = zeros(ntrials,1);
    for trId=1:ntrials
        cstart = TrialStart(trId);
        cstop = cstart + trial_dur - 1;
        TrialData(:,:,trId) = movavg_alpha(cstart:cstop,:);
        tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
    end
    
    %% Baseline extraction
    minFix_dur = min(FixStop - FixStart);
    Reference = NaN(minFix_dur, nchannels, ntrials);
    for trId=1:ntrials
        cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
        cstop = cstart+ minFix_dur - 1;
        Reference(:,:,trId) = movavg_alpha(cstart:cstop,:);
    end
    
    %% Compute ERD and LogBandPOwer
    Baseline = repmat(mean(Reference),[size(TrialData,1) 1 1]);
    %ERD = log(TrialData./Baseline);
    ERD = log(TrialData);      %Logband
    
    %%% Visualization %%%
    %% DIFFERENCE BETWEEN TASKS (IMAGESC)
    % img = (mean(ERD(:,:,tCk == classes(2)), 3) - mean(ERD(:,:,tCk == classes(1)), 3))';
    % img = img(find(~strcmp(channels_label, '')),:);
    % 
    % xcue = floor(mean(cue_pos - fix_pos));
    % xcf  = floor(mean(feedb_pos - fix_pos));
    % 
    %     figure();
    %     imagesc(img);
    %     colorbar;
    %     set(gca, 'CLim', [-0.6 0.6])
    %     line([xcue xcue], [0 40], 'Color', 'black', 'LineWidth', 2);
    %     line([xcf xcf], [0 40], 'Color', 'black', 'LineWidth', 2);
    %     xlabel('Time [s]')
    %     ylabel('Channels');
    %     xticks(0:events.SampleRate:minFix_dur);
    %     xticklabels((0:events.SampleRate:minFix_dur)/events.SampleRate);
    %     yticks(1:numel(channels_label(find(~strcmp(channels_label, '')))));
    %     yticklabels(channels_label(find(~strcmp(channels_label, ''))));
    
    %% TOPOPLOTS
    chanlocs_label = {chanlocs.labels};
    fixPeriod = [1/events.SampleRate 2]*events.SampleRate;
    cuePeriod = [2 3]*events.SampleRate;
    cfPeriod = [3 4; 4 5; 5 6; 6 (trial_dur/events.SampleRate)]*events.SampleRate;
    
    %Select channels
    recorded_channels = {'P3', 'Pz', 'P4', 'POz', 'O1', 'O2','P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'Oz'};
    %recorded_channels = {'Fp1','Fp2','F3','Fz','F4','FC1','FC2','C3','Cz','C4','CP1','CP2','P3','Pz','P4','POz','O1','O2','F1','F2','FC3','FCz','FC4','C1','C2','CP3','CPz','CP4','P5','P1','P2','P6','PO5','PO3','PO4','PO6','PO7','PO8','Oz'};
    % Trova gli indici dei canali registrati nella struct chanlocs.
    
    % figure()
    %% fix
    % dataFix_1 = mean(mean(ERD(fixPeriod(1):fixPeriod(2), :, tCk == classes(1)), 3), 1);
    % dataFix_2 = mean(mean(ERD(fixPeriod(1):fixPeriod(2), :, tCk == classes(2)), 3), 1);
    % dataFix = dataFix_2 - dataFix_1;
    % c_fix = zeros(64, 1);
    % for i=1:length(chanlocs_label)
    %     for j = 1:length(recorded_channels)
    %         if strcmpi(chanlocs_label{i}, recorded_channels{j})
    %             if ~isnan(dataFix(j))
    %                 c_fix(i) = dataFix(j);
    %             else
    %                 c_fix(i) = 0;
    %             end
    % 
    %         end
    %     end
    % end
    % 
    % subplot(2,3,1)
    % topoplot(squeeze(c_fix),chanlocs,'headrad', 'rim', 'maplimits', [-max(abs(c_fix)) max(abs(c_fix))]);
    % axis image;
    % title(['ERD/ERS fix (band [' num2str(band(1)) '-' num2str(band(2)) ']) - bottom right - bottom left']);
    % colorbar;
    
    %% cue
    % dataCue_1 = mean(mean(ERD(cuePeriod(1):cuePeriod(2), :, tCk == classes(1)), 3), 1);
    % dataCue_2 = mean(mean(ERD(cuePeriod(1):cuePeriod(2), :, tCk == classes(2)), 3), 1);
    % dataCue = dataCue_2 - dataCue_1;
    % c_cue = zeros(64, 1);
    % for i=1:length(chanlocs_label)
    %     for j = 1:length(recorded_channels)
    %         if strcmpi(chanlocs_label{i}, recorded_channels{j})
    %             if ~isnan(dataCue(j))
    %                 c_cue(i) = dataCue(j);
    %             else
    %                 c_cue(i) = 0;
    %             end
    % 
    %         end
    %     end
    % end
    % subplot(2,3,1)
    % topoplot(squeeze(c_cue),chanlocs,'headrad', 'rim', 'maplimits', [-max(abs(c_cue)) max(abs(c_cue))]);
    % axis image;
    % title(['ERD/ERS cue (band [' num2str(band(1)) '-' num2str(band(2)) ']) - bottom right - bottom left']);
    % colorbar;
    
    %% continuos feedback per intervals
    % for idx_cf = 1:size(cfPeriod,1)
    %     c_cfPeriod = cfPeriod(idx_cf,:);
    %     dataCf_1 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(1)), 3), 1);
    %     dataCf_2 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(2)), 3), 1);
    %     dataCf = dataCf_2 - dataCf_1;
    %     c_feedb = zeros(64, 1);
    %     for i=1:length(chanlocs_label)
    %         for j = 1:length(recorded_channels)
    %             if strcmpi(chanlocs_label{i}, recorded_channels{j})
    %                 if ~isnan(dataCf(j))
    %                 c_feedb(i) = dataCf(j);
    %                 else
    %                 c_feedb(i) = 0;
    %                 end
    % 
    %             end
    %         end
    %     end
    % 
    %     subplot(2,3,idx_cf+1)
    %     topoplot(squeeze(c_feedb),chanlocs,'headrad', 'rim', 'maplimits', [-max(abs(c_feedb)) max(abs(c_feedb))]);
    %     axis image;
    %     title(['ERD/ERS (band [' num2str(band(1)) '-' num2str(band(2)) ']) -- br - bl -- cf from ' num2str(ceil((cfPeriod(idx_cf,1) - cfPeriod(1,1))/events.SampleRate))...
    %         's to ' num2str(ceil((cfPeriod(idx_cf,2) - cfPeriod(1,1))/events.SampleRate)) 's']);
    %     colorbar;
    % end
    
    %% continuos feedback all
    c_cfPeriod = [3*events.SampleRate trial_dur];
    dataCf_1 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(1)), 3), 1);
    dataCf_2 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(2)), 3), 1);
    dataCf = dataCf_2 - dataCf_1;
      c_feedb = zeros(64, 1);
        for i=1:length(chanlocs_label)
            for j = 1:length(recorded_channels)
                if strcmpi(chanlocs_label{i}, recorded_channels{j})
                    if ~isnan(dataCf(j))
                    c_feedb(i) = dataCf(j);
                    else
                    c_feedb(i) = 0;
                    end
                end
            end
        end
    figure(1);
    subplot(2,3,f_idx)
    topoplot(squeeze(c_feedb), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_feedb)) max(abs(c_feedb))]);
    axis image;
    title(['ERD/ERS (band [' num2str(sel_band(1)) '-' num2str(sel_band(2)) ']) -- br - bl -- cf from 0' ...
         's to ' num2str(ceil((c_cfPeriod(2) - c_cfPeriod(1))/events.SampleRate)) 's']);
    colorbar;
    sgtitle(['Logband Continuos Feedback Subj: ' subject]);
    
    figure(2)
    whole_erd1 = mean(mean(ERD(:, :, tCk == classes(1)), 3), 1);
    whole_erd2 = mean(mean(ERD(:, :, tCk == classes(2)), 3), 1);
    whole_erd = whole_erd2-whole_erd1;
    total_logb = zeros(64, 1);
        for i=1:length(chanlocs_label)
            for j = 1:length(recorded_channels)
                if strcmpi(chanlocs_label{i}, recorded_channels{j})
                    if ~isnan(dataCf(j))
                    total_logb(i) = dataCf(j);
                    else
                    total_logb(i) = 0;
                    end
                end
            end
        end
    
    subplot(2,3,f_idx)
    topoplot(squeeze(total_logb), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(total_logb)) max(abs(total_logb))]);
    axis image;
    colorbar;
    sgtitle(['Logband Total Signal Subj: ' subject]);


    figure(3)
    cuetocf_erd1 = mean(mean(ERD(cuePeriod(1):cfPeriod(2), :, tCk == classes(1)), 3), 1);
    cuetocf_erd2 = mean(mean(ERD(cuePeriod(1):cfPeriod(2), :, tCk == classes(2)), 3), 1);
    cuetocf_erd = cuetocf_erd2-cuetocf_erd1;
    cuetoc_feed = zeros(64, 1);
        for i=1:length(chanlocs_label)
            for j = 1:length(recorded_channels)
                if strcmpi(chanlocs_label{i}, recorded_channels{j})
                    if ~isnan(dataCf(j))
                    cuetoc_feed(i) = dataCf(j);
                    else
                    cuetoc_feed(i) = 0;
                    end
                end
            end
        end
    subplot(2,3,f_idx)
    topoplot(squeeze(cuetoc_feed), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(cuetoc_feed)) max(abs(cuetoc_feed))]);
    axis image;
    colorbar;
    sgtitle(['Logband Cue to Cf Subj: ' subject]);


end
