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
subject = ['g2'];
lap_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Laplacian\lap_39ch_CVSA.mat';
chanlocs_path = 'C:\Users\User\Desktop\MATLAB\CVSA\Chanlocs\new_chanlocs64.mat';
path = ['C:\Users\User\Desktop\MATLAB\CVSA\records\' subject '\gdf'];
% ubuntu path = ['/home/riccardo/test_ws/records/' subject '/matselected_Trials];
classLb = {'Task 1','Task 2'};
classes = [730,731];

gdf_files = dir(fullfile(path, '*.gdf'));
load(lap_path);
load(chanlocs_path);

s=[]; events = struct('TYP',[],'POS',[],'SampleRate',512,'DUR',[]);
 for i=1:length(gdf_files)
        file = fullfile(path, gdf_files(i).name);
        [curr_s,h] = sload(file);
        curr_s = curr_s(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali
        slap = curr_s*lap;
        s = slap;
        curr_h = h.EVENT;
        events.TYP = curr_h.TYP;
        events.DUR = curr_h.DUR;
        events.POS = curr_h.POS;
    % concateno eventi
    % events.TYP = cat(1, events.TYP, curr_h.TYP);
    % events.DUR = cat(1, events.DUR, curr_h.DUR);
    % events.POS = cat(1, events.POS, curr_h.POS + size(s, 1));
    % s = cat(1, s, curr_s);
        
        %% Check valori del canale EOG
        % i dati sono nell'ordine dei microVolt: 2.02*10^4*10^-6
        targetElectrode = 'EOG';
        eog_ch= find(strcmp(channels_label, targetElectrode));
        idx_discard = find(abs(s(:,eog_ch))>2.5e+04); %threshold a 25mV

        %% Create Vector labels
        [nsamples,nchannels] = size(s);
        [feedb_pos, feedb_dur, fix_dur, fix_pos, cue_dur, cue_pos, ntrials] = extract_info_label(events, 781, 786, [730 731]);

        %% Extract trial data
        [TrialStart, TrialStop, FixStart, FixStop, Ck, Tk] = extract_trial_info(s, events, fix_pos, fix_dur, feedb_pos, feedb_dur, cue_pos, ntrials);

        %Se i valori di s nel canale EOG sono <25mV tebgo il trial e lo
        %concateno con il resto del segnale, altrimenti lo scarto e passo a
        %quello dopo o li pongo a zero.
        %Posso usare il find per trovare gli indici
        
        % Con g2 tutte le run sono oltre il threshold, negli altri apposto
        if isempty(idx_discard)
            disp('No trial to be discarded')
        else
              tk_discard = Tk(idx_discard);
              trial_to_discard = unique(nonzeros(tk_discard));
              disp(['Trial to de discarded for calibration run ', num2str(i)])
              disp(trial_to_discard);

        end



 end