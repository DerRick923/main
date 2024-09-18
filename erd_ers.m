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
load(chanlocs_path);
features_file = ['/home/riccardo/test_ws/records/' c_subject '/dataset/selected_features.mat'];
features = load(features_file);
%selbands = features.selectedFeatures(:,2);
selchs = features.selectedFeatures(:,1);
channles_string = string(channels_label);

classLb = {'CVSA Left','CVSA Right'};
classes = [730,731];
nclasses = length(classes);



%files = dir(fullfile(path, '*.mat'));
files = dir(fullfile(path, '*.gdf'));

band = {[8 10], [10 12], [12 14], [14 16], [16 18]};
nbands = length(band);

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]); Rk=[];
for i=1:length(files)
    file = fullfile(path, files(i).name);

    %load(file);
    [signal,header] = sload(file);
    
    curr_s = signal(:,1:39);
    %slap = curr_s*lap;
    %curr_s = slap;
    curr_h = header.EVENT;
    
    if strcmp(test_typ, "calibration")
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
new_Rk = []; new_Ck = [];
DataperTrial = NaN(trial_dur,nchannels,nbands,ntrials);
tCk = zeros(ntrials,1);
for trId=1:ntrials
    cstart = TrialStart(trId);
    cstop = cstart + trial_dur - 1;
    DataperTrial(:,:,:,trId) = s_processed(cstart:cstop,:,:);

    c_Rk = Rk(cstart:cstop,1);
    new_Rk = cat(1,new_Rk,c_Rk);
    c_Ck = Ck(cstart:cstop);
    new_Ck = cat(1,new_Ck,c_Ck);
    tCk(trId) = unique(nonzeros(Ck(cstart:cstop)));
end

%% Baseline extraction for each trial
minFix_dur = min(FixStop - FixStart);
Reference = NaN(minFix_dur, nchannels, nbands, ntrials);
for trId=1:ntrials
    cstart = FixStart(trId); %=TrialStart(trId) o fix_pos(trId)
    cstop = cstart+ minFix_dur - 1;
    Reference(:,:,:,trId) = s_processed(cstart:cstop,:,:);
end
Baseline = repmat(mean(Reference),[size(DataperTrial,1) 1 1 1]);

%% Compute ERD and LogBandPOwer [samples x channels x bands] con tutti i trial
ERD = log(DataperTrial./Baseline);

%Segna canali di interesse e plotta ERD e ERS
%Quali tempi?
t = linspace(0,events.SampleRate,trial_dur);
for chan=1:length(selchs)
    chan_idx = find(selchs{chan}==channles_string);
    figure(chan)
    for i=1:nbands
            c_band = band{i};
            subplot(2,3,i)
            cERD_1 = mean(ERD(:,chan_idx,i,tCk==classes(1)),4);
            cERD_2 = mean(ERD(:,chan_idx,i,tCk==classes(2)),4);
            plot(t,cERD_1)
            grid on
            hold on
            plot(t,cERD_2)
            hold off
            axis tight
            legend('Bottom left (730)','Bottom right (731)')
            title(['[' num2str(c_band(1)) '-' num2str(c_band(2)) ']'])
            sgtitle(['ERD/ERS Subj: ' c_subject ' Channel: ' selchs{chan}])

    
    end
end