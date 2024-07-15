close all
clear all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\eeglab2024.0'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))


channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
       '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

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

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]);

band = {[6 9], [9 12], [8 14], [12 15], [15 18], [18 21]};

for f_idx=1:length(band)
    sel_band = band{f_idx}; %[hz]
    for runId=1:length(matfiles)
        file = fullfile(path, matfiles(runId).name);
        load(file);
        curr_s = signal(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali

        %slap = curr_s*lap;
        %s = slap;
        s = curr_s;
        curr_h = header.EVENT;
        events.TYP = curr_h.TYP;
        events.DUR = curr_h.DUR;
        events.POS = curr_h.POS;
    % concateno eventi
    % events.TYP = cat(1, events.TYP, curr_h.TYP);
    % events.DUR = cat(1, events.DUR, curr_h.DUR);
    % events.POS = cat(1, events.POS, curr_h.POS + size(s, 1));
    % s = cat(1, s, curr_s);

        %% Create Vector labels
        [nsamples,nchannels] = size(s);
        [feedb_pos, feedb_dur, fix_dur, fix_pos, cue_dur, cue_pos, ntrials] = extract_info_label(events, 781, 786, [730 731]);

        %% Extract trial data
        [TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(s, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, ntrials);

        %% Processing offline
        t_window = 1; %[s]
        windowSize = events.SampleRate*t_window;
        filtOrder = 4;
        signal_processed = data_processing(s,nchannels,events.SampleRate,sel_band,filtOrder,t_window);


        %% Trial extraction
        trial_dur = min(TrialStop-TrialStart);
        TrialData = NaN(trial_dur,nchannels,ntrials);
        tCk = zeros(ntrials,1);
        for trId=1:ntrials
            cstart = TrialStart(trId);
            cstop = cstart + trial_dur - 1;
            TrialData(:,:,trId) = signal_processed(cstart:cstop,:);
            tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
        end

        %% Baseline extraction
        minFix_dur = min(FixStop - FixStart);
        Reference = NaN(minFix_dur, nchannels, ntrials);
        for trId=1:ntrials
            cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
            cstop = cstart+ minFix_dur - 1;
            Reference(:,:,trId) = signal_processed(cstart:cstop,:);
        end

        %% Compute ERD and LogBandPOwer
        Baseline = repmat(mean(Reference),[size(TrialData,1) 1 1]);
        %ERD = log(TrialData./Baseline);
        ERD = log(TrialData);      %Logband


        %% VISUALIZATION
        fixPeriod = [1/events.SampleRate 2]*events.SampleRate;
        cuePeriod = [2 3]*events.SampleRate;
        cfPeriod = [3 4; 4 5; 5 6; 6 (trial_dur/events.SampleRate)]*events.SampleRate;
        chanlocs_label = {chanlocs.labels};

        %Select channels
        recorded_channels = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
        %recorded_channels = {'Fp1','Fp2','F3','Fz','F4','FC1','FC2','C3','Cz','C4','CP1','CP2','P3','Pz','P4','POz','O1','O2','F1','F2','FC3','FCz','FC4','C1','C2','CP3','CPz','CP4','P5','P1','P2','P6','PO5','PO3','PO4','PO6','PO7','PO8','Oz'};
% Trova gli indici dei canali registrati nella struct chanlocs.
% recorded_indices = zeros(1,length(recorded_channels));
% for i = 1:length(recorded_channels)
%     recorded_indices(i) = find(strcmp({chanlocs.labels},recorded_channels{i}));
% end

        % TOPOPLOT
        % Per ogni banda di frequenze per ogni run
        % All cf for the two tasks
        figure(f_idx)
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

        % Computation for the different tasks: c_feedb_1 = zeros(64, 1);
        %     for idx=1:length(chanlocs_label)
        %         for j = 1:length(recorded_channels)
        %             if strcmpi(chanlocs_label{idx}, recorded_channels{j})
        %                 if ~isnan(dataCf_1(j))
        %                 c_feedb_1(idx) = dataCf_1(j);
        %                 else
        %                 c_feedb_1(idx) = 0;
        %                 end
        % 
        %             end
        %         end
        %     end
        % 
        % c_feedb_2 = zeros(64, 1);
        %     for idx=1:length(chanlocs_label)
        %         for j = 1:length(recorded_channels)
        %             if strcmpi(chanlocs_label{idx}, recorded_channels{j})
        %                 if ~isnan(dataCf_2(j))
        %                 c_feedb_2(idx) = dataCf_2(j);
        %                 else
        %                 c_feedb_2(idx) = 0;
        %                 end
        % 
        %             end
        %         end
        %     end
        %     subplot(1,2,1)
        %     topoplot(squeeze(c_feedb_1), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(dataCf_1)) max(abs(dataCf_1))]);
        %     axis image;
        %     % title(['ERD/ERS (band [' num2str(band(1)) '-' num2str(band(2)) ']) -- br - bl -- cf from 0' ...
        % %     's to ' num2str(ceil((c_cfPeriod(2) - c_cfPeriod(1))/sampleRate)) 's']);
        %     colorbar;
        %     subplot(1,2,2)
        %     topoplot(squeeze(c_feedb_2), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_feedb_2)) max(abs(c_feedb_2))]);
        %     axis image;
        %     % title(['ERD/ERS (band [' num2str(band(1)) '-' num2str(band(2)) ']) -- br - bl -- cf from 0' ...
        % %     's to ' num2str(ceil((c_cfPeriod(2) - c_cfPeriod(1))/sampleRate)) 's']);
        %     colorbar;

        subplot(2,ceil(length(matfiles)/2),runId)
        topoplot(squeeze(c_feedb), chanlocs,'electrodes','labelpoint','headrad', 'rim', 'maplimits', [-max(abs(c_feedb)) max(abs(c_feedb))]);
        axis image;
        title(['Run ' num2str(runId) ', Band [' num2str(sel_band(1)) '-' num2str(sel_band(2)) ']']);
        colorbar;
        sgtitle(['ERD Continuos Feedback in [' num2str(sel_band(1)) '-' num2str(sel_band(2)) ']'])


        %% cue+cf
        figure(length(band)+f_idx)
        dataCue_1 = mean(mean(ERD(cuePeriod(1):cfPeriod(2), :, tCk == classes(1)), 3), 1);
        dataCue_2 = mean(mean(ERD(cuePeriod(1):cfPeriod(2), :, tCk == classes(2)), 3), 1);
        dataCue = dataCue_2 - dataCue_1;
        c_cue = zeros(64, 1);
        for i=1:length(chanlocs_label)
            for j = 1:length(recorded_channels)
                if strcmpi(chanlocs_label{i}, recorded_channels{j})
                    if ~isnan(dataCue(j))
                        c_cue(i) = dataCue(j);
                    else
                        c_cue(i) = 0;
                    end
                end
            end
        end
        subplot(2,ceil(length(matfiles)/2),runId)
        topoplot(squeeze(c_cue), chanlocs,'electrodes','labelpoint','headrad', 'rim', 'maplimits', [-max(abs(c_feedb)) max(abs(c_feedb))]);
        axis image;
        title(['Run ' num2str(runId) ', Band [' num2str(sel_band(1)) '-' num2str(sel_band(2)) ']']);
        colorbar;
        sgtitle(['ERD from Cue to CF in [' num2str(sel_band(1)) '-' num2str(sel_band(2)) ']'])

    end
end
