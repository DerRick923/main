% compute fischer score matrix with dimensions [freq_bands x channels x runs]
close all
clear all
clc

addpath(genpath('/home/riccardo/Desktop/CVSA'))
addpath(genpath('/home/riccardo/lib/cnbi-smrtrain'))
addpath(genpath('/home/riccardo/test_ws'))
addpath(genpath('/home/riccardo/Desktop/eeglab2024.0'))

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
       '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};
 
%all_channels = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%       'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};


% file info
c_subject = 'h7';
prompt = 'Enter "calibration" or "evaluation": ';
test_typ = input(prompt, 's');

%path = ['/home/riccardo/test_ws/records/' c_subject '/matselected_Trials];
path = ['/home/riccardo/test_ws/records/' c_subject '/gdf/' test_typ];
%path = ['/home/riccardo/test_ws/records/' c_subject '/gdf/calibration'];
chanlocs_path = '/home/riccardo/Desktop/CVSA/Chanlocs/new_chanlocs64.mat';

classLb = {'Bottom right','Bottom left'};
classes = [730,731];
nclasses = length(classes);
bands = {[8 10], [10 12], [12 14], [14 16], [16 18], [8 14]};
nbands = length(bands);
sampleRate = 512;

load(chanlocs_path);
chanlocs_label = {chanlocs.labels};
%files = dir(fullfile(path, '*.mat'));
files = dir(fullfile(path, '*.gdf'));
for idx_f=1:length(files)
    file = fullfile(path, files(idx_f).name);

    %load(file);
    [signal,header] = sload(file);
    
    s = signal(:,1:39);
    %slap = curr_s*lap;
    %curr_s = slap;
    events = header.EVENT;
    
    if strcmp(test_typ, "calibration")
        start = find(events.TYP == 1,1,'first');
        events.TYP = events.TYP(start:end);
        events.POS = events.POS(start:end);
        events.DUR = events.DUR(start:end);
    end

    % Create Vector labels
    [nsamples,nchannels] = size(s);
    [feedb_pos, feedb_dur, fix_pos, fix_dur, cue_pos, cue_dur, ntrials] = extract_info_label(events, 781, 786, [730 731]);

    % Extract trial data
    [TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(s, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, ntrials);
    
    for idx_band=1:nbands
        band = bands{idx_band};
        % Processing offline
        t_window = 1; %[s]
        windowSize = events.SampleRate*t_window;
        filtOrder = 4;
        signal_processed = data_processing(s,nchannels,events.SampleRate,band,filtOrder,t_window);
        
        % Trial extraction
        trial_dur = min(TrialStop-TrialStart);
        TrialData = NaN(trial_dur,nchannels,ntrials);
        tCk = zeros(ntrials,1);
        for trId=1:ntrials
            cstart = TrialStart(trId);
            cstop = cstart + trial_dur - 1;
            TrialData(:,:,trId) = signal_processed(cstart:cstop,:);
            tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
        end

        % Baseline extraction
        minFix_dur = min(FixStop - FixStart);
        Reference = NaN(minFix_dur, nchannels, ntrials);
        for trId=1:ntrials
            cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
            cstop = cstart+ minFix_dur - 1;
            Reference(:,:,trId) = signal_processed(cstart:cstop,:);
        end

        % Compute ERD and LogBandPOwer
        Baseline = repmat(mean(Reference),[size(TrialData,1) 1 1]);
        %ERD = log(TrialData./Baseline);
        ERD = log(TrialData);      %Logband

        % Visualization
        %Scannerizzo il trial ogni mezzo secondo per vedere l'evoluzione
        %delle logband nelle varie fasi
        period1 = 1:sampleRate/2:trial_dur;
        period2 = period1(2):sampleRate/2:trial_dur;
        period2 = cat(2, period2, trial_dur);
        period = cat(1, period1, period2);
        
        % show for each period the topoplot
        figure();
        for idx_period=1:size(period, 2)
            data_1 = mean(mean(ERD(period(1, idx_period):period(2, idx_period), :, tCk == classes(1)), 3), 1);
            data_2 = mean(mean(ERD(period(1, idx_period):period(2, idx_period), :, tCk == classes(2)), 3), 1);
            data = data_2 - data_1;
            chanlocs_data = zeros(size(chanlocs_label,2), 1);
            for i=1:length(chanlocs_label)
                for j = 1:nchannels
                    if strcmpi(chanlocs_label{i}, channels_label{j})
                        if ~isnan(data(j))
                            chanlocs_data(i) = data(j);
                        else
                            chanlocs_data(i) = 0;
                        end
        
                    end
                end
            end

            subplot(2, ceil((size(period, 2) + 1)/2), idx_period);
            topoplot(squeeze(chanlocs_data), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(chanlocs_data)) max(abs(chanlocs_data))]);
            axis image;
            colorbar;
            if period(2, idx_period) - 1 <= sampleRate*2
                title(['fixation: ' num2str(period(1, idx_period)/sampleRate) '-' num2str(period(2, idx_period)/sampleRate) 's']);
            elseif period(2, idx_period) - 1 <= sampleRate*3
                title(['cue: ' num2str(period(1, idx_period)/sampleRate) '-' num2str(period(2, idx_period)/sampleRate) 's']);
            else
                title(['cf: ' num2str(period(1, idx_period)/sampleRate) '-' num2str(period(2, idx_period)/sampleRate) 's']);
            end
        end
    
        % show the topoplot for all the cf
        c_cfPeriod = [3*sampleRate trial_dur];
        dataCf_1 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(1)), 3), 1);
        dataCf_2 = mean(mean(ERD(c_cfPeriod(1):c_cfPeriod(2), :, tCk == classes(2)), 3), 1);
        dataCf = dataCf_2 - dataCf_1;
        c_cf = zeros(64,1);
        for i=1:length(chanlocs_label)
            for j = 1:nchannels
                if strcmpi(chanlocs_label{i}, channels_label{j})
                    if ~isnan(dataCf(j))
                        c_cf(i) = dataCf(j);
                    else
                        c_cf(i) = 0;
                    end
                end
            end
        end
        
        subplot(2, ceil((size(period, 2) + 1)/2), idx_period + 1);
        topoplot(squeeze(c_cf), chanlocs, 'headrad', 'rim', 'maplimits', [-max(abs(c_cf)) max(abs(c_cf))], 'electrodes', 'labelpoint');
        axis image;
        title(['all the cf: ' num2str(c_cfPeriod(1)/sampleRate) '-' num2str(c_cfPeriod(2)/sampleRate) 's']);
        colorbar;
        
        all_title = [files(idx_f).name '| br-bl | band: [' num2str(band(1)) ',' num2str(band(2)) ']'];
        
        sgtitle(all_title)
    end
end