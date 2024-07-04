clear all
close all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))

channels_label = {'', '', '', '', '', '', '', '', '', '', '', '', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', '', ...
        '', '', '', '', '', '', '', '', '', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% channels_label = {'FP1', 'FP2', 'F3', 'FZ', 'F4', 'FC1', 'FC2', 'C3', 'CZ', 'C4', 'CP1', 'CP2', 'P3', 'PZ', 'P4', 'POZ', 'O1', 'O2', 'EOG', ...
%         'F1', 'F2', 'FC3', 'FCZ', 'FC4', 'C1', 'C2', 'CP3', 'CP4', 'P5', 'P1', 'P2', 'P6', 'PO5', 'PO3', 'PO4', 'PO6', 'PO7', 'PO8', 'OZ'};

% file info
subject = 'c7';
date = '240503';
lap_path = "C:\Users\User\Desktop\MATLAB\CVSA\Laplacian\lap_39ch_CVSA.mat";
chanlocs_path = "C:\Users\User\Desktop\MATLAB\CVSA\Chanlocs\new_chanlocs64.mat";
path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' subject '\mat_selectedTrials'];
% path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\ ' subject '\gdf'];
% ubuntu path = ['/home/riccardo/test_ws/records/ ' subject '/matselected_Trials];

classes = [730, 731];
band = [8 14];
samplerate = 512;
nchannels = length(channels_label);

matfiles = dir(fullfile(path,'*.mat'));
load(lap_path);

%% CONCATENAZIONE
events = struct('POS', [], 'DUR', [], 'TYP', []);
total_signal = [];

for i=1:length(matfiles)
    file = fullfile(path, matfiles(i).name);
    load(file);
    curr_s = signal(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali
    curr_h = header.EVENT;
    slap = curr_s*lap;

    total_signal = cat(1,total_signal, signal);
    events.POS = cat(1,events.POS,curr_h.POS + length(total_signal));
    events.DUR = cat(1,events.DUR,curr_h.DUR);
    events.TYP = cat(1,events.TYP,curr_h.TYP);

    info.classes = classes;
    info.sampleRate = samplerate;
    info.lap = lap;
    info.selChan = [];
    info.selFreq = [];

end
