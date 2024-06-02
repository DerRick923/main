clear all
close all
clc

addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\matlab'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t200_FileAccess'))
addpath(genpath('C:\Users\User\Desktop\biosig-2.5.1-Windows-64bit\biosig-2.5.1-Windows-64bit\share\matlab\t250_ArtifactPreProcessingQualityControl'))
addpath(genpath('C:\Users\User\Desktop\MATLAB\CVSA'))

matfiles = ["c7.20240417.111204.calibration.cvsa_lbrb.mat","c7.20240417.111557.calibration.cvsa_lbrb.mat"];
load('lap_39ch_CVSA.mat')

filename1 = 'c7.20240417.111204.calibration.cvsa_lbrb.mat';
filename2 = 'c7.20240417.111557.calibration.cvsa_lbrb.mat';
files = {filename1,filename2};

for i=1:length(matfiles)
    cdata = load(matfiles(i));
    cfilename = files{i};
    curr_s = cdata.signal(:,1:39);    %ho 39 canali e la matrice ha 40 colonne, quindi seleziono solo le colonne riferite ai canali
    curr_h = cdata.header;
    slap = curr_s*lap;
   
    %% Calcolo PSD (da fare per ogni segnale)
    wlength = 0.5; % seconds. Length of the external window
    pshift = 0.25; % seconds. Shift of the internal windows
    wshift = 0.0625; % seconds. Shift of the external window
    samplerate = 512;   % curr_h.SampleRate (??)
    mlength = 1; % seconds
    [psd, f] = proc_spectrogram(slap, wlength, wshift, pshift, samplerate, mlength);

    %Recompute the EVENT POS and DUR with respect to PSD windows
    winconv = 'backward';
    POS = proc_pos2win(curr_h.POS', wshift*samplerate, winconv, wlength*samplerate);
    DUR = proc_pos2win(curr_h.DUR', wshift*samplerate, winconv, wlength*samplerate);
    event = struct('TYP',curr_h.TYP','POS',POS,'SampleRate',samplerate,'DUR',DUR);

% Per salvare i dati di ogni file e riutilizzarli in altri script.
% data = struct();
% data.PSD = psd;
% data.selFrequencies = f;
% data.event = event;

[~, pfilename] = fileparts(cfilename);
sfilename = [pfilename '_filt.mat'];
save(sfilename,'psd','f','event','samplerate');

end

%togliere kaplaciano e provare con filtraggio e moving average
%calolare band power al post di erd


