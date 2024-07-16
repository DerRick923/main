function result = eye_movement_check(data,channels_label,threshold)
    %Fai il check di Fp1,Fp2,EOG per valori > 25mv
    % Se sono maggiori il trial è da scartare, altrimenti tieni il trial e vai
    % avanti. I dati sono nell'ordine dei microVolt: 2.02*10^4*10^-6

    % INPUT:
    %   data: eeg signal per trial
    %   channels_label: cell vector with all the channels used for
    %   recordings
    %   threshold: value to define if the trial should be discarded

    % OUTPUT:
    %   result: returns true if all values for the specific channels are > threshold, otherwise return false
    result = false;
    target_electrodes = {'FP1','FP2','EOG'};
    eog_ch= find(ismember(channels_label, target_electrodes));
    
        for i=1:length(eog_ch)
            if any(abs(data(:,eog_ch(i)))>threshold) %threshold a 25mV
                result = true;
                return;
            end
        end       
end
