clear all
close all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))

% channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
%        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% file info
subject = 'c7';
date = '240503';
lap_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Laplacian\lap_39ch_CVSA.mat';
chanlocs_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Chanlocs\new_chanlocs64.mat';
path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' subject '\mat_selectedTrials'];
% path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\ ' subject '\gdf'];
% ubuntu path = ['/home/riccardo/test_ws/records/ ' subject '/matselected_Trials];

classes = [730, 731];
band = [8 14];
nchannels = length(channels_label);

matfiles = dir(fullfile(path, '*.mat'));
load(lap_path);

%% INFO
wlength = 0.5; % seconds. Length of the external window
pshift = 0.25; % seconds. Shift of the internal windows
wshift = 0.0625; % seconds. Shift of the external window
samplerate = 512;   % curr_h.SampleRate (??)
mlength = 1; % seconds

%% CONCATENAZIONE
PSD = []; events = struct('POS', [], 'DUR', [], 'TYP', [], 'info',[]);
for i=1:length(matfiles)
    file = fullfile(path, matfiles(i).name);
    load(file);
    signal = signal(:,1:nchannels);
    curr_h = header.EVENT;
    
    %% compute psd for each file/run
    signal = signal * lap;
    [psd,f] = proc_spectrogram(signal,wlength,wshift,pshift,samplerate,mlength);
    psd = log(psd);
    winconv = 'backward';
    curr_h.POS = proc_pos2win(curr_h.POS, wshift*samplerate, winconv, wlength*samplerate);
    curr_h.DUR = proc_pos2win(curr_h.DUR, wshift*samplerate, winconv, wlength*samplerate);

    psd = squeeze(mean(psd(:, find(f >= band(1) & f <= band(2)),:),2));
    
    % Concatenate every psd
    events.POS = cat(1,events.POS,curr_h.POS + size(PSD,1));
    events.DUR = cat(1,events.DUR,curr_h.DUR);
    events.TYP = cat(1,events.TYP,curr_h.TYP);
    PSD = cat(1,PSD,psd);

    info.wlength = wlength;
    info.pshift = pshift;
    info.wshift = wshift;
    info.samplerate = samplerate;
    info.mlength = mlength;
end

%% Extract info data
[feedb_pos, feedb_dur, fix_dur, fix_pos, cue_dur, cue_pos, trials] = extract_info_label(events, 781, 786, [730 731]);

%% Extract trial data
[TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(PSD, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, trials);

nwindows = size(PSD,1);
%nfreq = size(PSD,2); con la media delle frequenze non c'è più nfreq
nchannels = size(PSD,2);

% Create Activity matrix
trial_dur = min(TrialStop-TrialStart);
Activity = NaN(trial_dur,nchannels,trials);
%Activity = NaN(trial_dur,nfreq,nchannels,trials);
tCk = zeros(trials,1);
for trId=1:trials
    cstart = fix_pos(trId);
    cstop = cstart + trial_dur - 1;
    Activity(:,:,trId) = PSD(cstart:cstop,:);
    tCk(trId) = unique(Ck(cstart:cstop));
end

% Compute the difference in the psd between the two tasks
  psd_1 = mean(Activity(:,:,tCk == classes(1)), 3);
  psd_2 = mean(Activity(:,:,tCk == classes(2)), 3);
  psd_diff = (psd_2 - psd_1);
  psd_diff = psd_diff'; 
   
  xcue = mean(cue_pos - fix_pos) +1;
  xcf = mean(feedb_pos - fix_pos) +1;

    figure();
    imagesc(psd_diff);
    line([xcue xcue], [1 39], 'Color', 'black', 'LineWidth', 2);
    line([xcf xcf], [1 39], 'Color', 'black', 'LineWidth', 2);
    colorbar;
    set(gca, 'CLim', [-0.6 0.6])
    xlabel('Time [s]')
    %xticks(1:size(psd_diff,2));
    %xticklabels(0:1/16:size(psd_diff,2)/16);
    ylabel('Channels');
    yticks(1:numel(channels_label));
    yticklabels(channels_label);
    title('Feature Map: psd');



