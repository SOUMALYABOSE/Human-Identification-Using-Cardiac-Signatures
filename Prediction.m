%% Important Copyright Information with regards to this scirpt:
% Part(s) of the Script has been taken from 2020 Script of: Sven Schellenberger, Kilin Shi (Copyright (C) 2020  Sven Schellenberger, Kilin Shi)
% Rest has been added under the CME-Master Thesis titled "Human Identification With Gender Knowledge Analysis From Cardiac Signatures and Health-related Information Extraction From Radar Based Biomarkers" 
% Conducted @ Home Automation Lab, Institute FAPS, FAU Erlangen-Nuremberg by Mr. Soumalya Bose
% Advisors: Mr. Jochen Bauer (Institute FAPS, FAU), Dr. -med. Tobias Steigleder (PallMeT, Univ. Hospital Erlangen, FAU), Professor Dr.-Ing. Jörg Franke (Institute FAPS, FAU), Professor Dr.-Ing. Georg Fischer (Institute for Electronics Engineering, FAU)
% Copyright (C) 2023-24  Mr. Soumalya Bose
%% Plot data
% This script generates exemplary plots of selected subject-ID(s) and
% scenario(s)
% Data is loaded, processed and then plotted

%% Clearing of Workspace
clear;
clc;

%% Init
addpath(genpath('utils'))

%% Preparation of Datasets

% GDN00XX - simply choose numbers or ranges from 01 to 30
  IDrange = [5 8 12 16 17 19 23 28 30]; % (Uncomment/Comment When Needed)
% IDrange = 8; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
% scenarios = {'Resting', 'Valsalva', 'Apnea', 'TiltUp', 'TiltDown'}; % (Uncomment/Comment When Needed)
  scenarios = {'Resting'}; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
 
%%%% Declaring Confusion Matrix Parameters
TP = 0;
TN = 0;
FP = 0;
FN = 0;

%%%% Set path to datasets 
% Datasets can be found on figshare
path = 'datasets'; % In this case "datasets" folder is in this scripts folder  

scrsz = get(groot,'ScreenSize'); % For plotting
%% Extract and plot the data

% Manually selected ECG channel for each subject (Position in ECG_CHANNEL vector corresponds to ID)
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];

for indx = 1:length(IDrange)
    %% Iterate all subject IDs
    ID = sprintf('GDN%04d',IDrange(indx));
    % fprintf('----------- Loading %s ------------\n', ID);

    for sz = 1:length(scenarios)
        %% Iterate all existing subject scenarios
        scenario = scenarios{sz};
        % fprintf('---- Scenario %s\n', scenario);
        
        % search file
        path_id = [path,'\',ID];
        files_synced_mat = dir([path_id,'\*.mat']);
        found = [];
        for j = 1:length(files_synced_mat)
            found = strfind(files_synced_mat(j).name,scenario);
            if ~isempty(found)
                % load file
                load([path_id,'\',files_synced_mat(j).name]);
                break;
            end
        end
        if isempty(found)
            fprintf('---- skipped\n');
            ref_autocorr_array = 0; % Customized by Soumalya Bose to avoid Technical Glitches
            ref_PS_w_norm = 0; % Customized by Soumalya Bose to avoid Technical Glitches
            continue
        end
        
        %% Prepare data
        output.(ID).(scenario) = struct;
        [~,~,~,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},0); % Ellipse fitting, compensation and distance reconstruction
        % Usage of elreko
        % [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1}(timestamp of dataset),0(Plot flag -> 1: plots on, 0: plots off));
        
        % Radar
        [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);
               
        %% Heartbeat detection
        radar_hs_states = getHsmmStates(radar_heartsound, fs_radar);
        radar_hs_locsR = statesToLocs(radar_hs_states, 1);      
        
        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar;        

        %% Calculate autocorrelation
        radar_mod_heartsound = radar_heartsound.*10^6;
        lags = -length(radar_mod_heartsound) + 1 : length(radar_mod_heartsound) - 1;
        start_index = 1;  % Specify the starting index for autocorrelation calculation
        end_index = min(start_index + 10000, length(radar_mod_heartsound));  % Specify the ending index, limit to a reasonable size
        [autocorr_values, lag_indices] = xcorr(radar_mod_heartsound(start_index:end_index), 'biased');
        ref_autocorr_array = autocorr_values;

        %% Calculate Power Spectrum
        
        % Pre-processing of Radar Data/Rel. Distance/Chest Displacement Signal      
        radar_dist_scaled = radar_dist*1000; % Scaling of Raw Radar Displacement Signal
        % Kaiser Window with β = 5.1 (0-4: Wide ML-L. Suppre SL-High Freq Resol. 4-10: Reverse+Less Spectral Leak/Artefact)
        w = kaiser(length(radar_dist_scaled), 5.1); 
        radar_dist_windowed = radar_dist_scaled.*w;
        Yw = fft(radar_dist_windowed);
        L = length(radar_dist_windowed);

        % Fequency Domain Transformation of Processed Radar Data/Rel. Distance/Chest Displacement Signal 
        PS2w = abs(Yw/L); % Double sampling plot
        PS1w = PS2w(1:floor(L/2)+1); % Single sampling plot. Used Floor func. as L is from length func. which returns double and not integer.
        PS1w(2:end-1) = 2*PS1w(2:end-1);
        f = fs_radar*(0:(L/2))/L;  % Frequency bins

        % Normalization of Magnitude Spectrum from 0.1 Hz to 1.0 Hz
        PS1_vals = interp1(f,PS1w,0.1:0.00004:1);
        max_PS1w = max(PS1_vals);
        PS_w_norm = PS1w/max_PS1w;
        ref_PS_w_norm = interp1(f,PS_w_norm,0.1:0.00004:0.8);

        %% Calculate Fundamental Frequency

        ref_fundamentalFrequency = interp1(PS_w_norm, f, 1);
        % fprintf("****Fundamental Frequency Position is: %f Hz. \n", ref_fundamentalFrequency);

        %% Calculate Second Breathing Harmonic

        start_pos = 2 * round(ref_fundamentalFrequency * 10) / 10;
        end_pos = start_pos + 0.1;
        breathHarmonic_freqRange_idx = (f >= start_pos) & (f <= end_pos);
        breathHarmonic_values = interp1(f, PS_w_norm, f(breathHarmonic_freqRange_idx));
        [second_breathHarmonic, index] = max(breathHarmonic_values);
        second_breathHarmonicFrequency = f(breathHarmonic_freqRange_idx);
        ref_second_breathHarmonicFrequency = second_breathHarmonicFrequency(index);
        % fprintf("****2nd Breathing Harmonic Position is: %f Hz. \n", ref_second_breathHarmonicFrequency);

        %% Calculate First Inter-Harmonic Distance

        ref_first_interHarmonicDist = ref_second_breathHarmonicFrequency - ref_fundamentalFrequency;
        % fprintf("****1st Inter-Harmonic Distance is: %f Hz. \n", ref_first_interHarmonicDist);

        %% Prediction

        % Calculating KL Divergence between Different Autocorrelation Functions cases
        if (scenario == "Resting")
            current_IDrange = [5 8 12 16 17 19 23 28 30]; % (Uncomment/Comment When Needed)
          % current_IDrange = [8 9]; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
            current_scenarios = {'Resting', 'Valsalva', 'Apnea'};

            for current_indx = 1:length(current_IDrange)
                current_ID = current_IDrange(current_indx);
                current_ID_full = sprintf('GDN%04d',current_ID);
                for current_sz = 1:length(current_scenarios)
                    current_scenario = current_scenarios{current_sz};
                    
                    %% Function Call
                    % Power Spectrum Function Call
                    [current_f, ff_PS_w_norm, current_PS_w_norm] = Ps_w_norm_calc(current_ID, current_scenario);

                    % Autocorrelation Function Call
                    current_autocorr_array = Autocorrl_calc(current_ID, current_scenario);

                        %% Power Spectrum Section                  
                        % Power Spectrum Array Padding
                        sizeDiff = abs(numel(ref_PS_w_norm) - numel(current_PS_w_norm)); 
                        if (numel(ref_PS_w_norm) < numel(current_PS_w_norm))
                            ref_PS_w_norm_padded = [ref_PS_w_norm zeros(1, sizeDiff)]; % Zero Padding: Essential for KL Divergence
                            % Conversion To Probability Distribution: Essential for KL Divergence
                            prob_ref_PS_w_norm = ref_PS_w_norm_padded ./ sum(ref_PS_w_norm_padded);
                            prob_current_PS_w_norm = current_PS_w_norm ./ sum(current_PS_w_norm); 
                        else
                            current_PS_w_norm_padded = [current_PS_w_norm zeros(1, sizeDiff)]; % Zero Padding: Essential for KL Divergence
                            % Conversion To Probability Distribution: Essential for KL Divergence
                            prob_ref_PS_w_norm = ref_PS_w_norm ./ sum(ref_PS_w_norm);
                            prob_current_PS_w_norm = current_PS_w_norm_padded ./ sum(current_PS_w_norm_padded);
                        end

                        % Compute the KL divergence of Power Spectrum
                        kl_divergence_ps = sum(prob_ref_PS_w_norm .* log2(prob_ref_PS_w_norm ./ prob_current_PS_w_norm), 'all');
                        normalized_divergence_ps = 10000*(kl_divergence_ps / (numel(prob_ref_PS_w_norm) + numel(prob_current_PS_w_norm)));
                        rounded_divergence_ps = round(normalized_divergence_ps, 2);  % Round to 2 decimal places
                        % fprintf("****Normalized Power Spectrum KL Divergence Score from Subject ID: %s Scenario: %s to Subject ID: %s Scenario: %s is: %f. \n", ID, scenario, current_ID_full, current_scenario, normalized_divergence_ps);
                        
                        %% Autocorrelation Section
                        % Autocorrelation Array Padding
                        sizeDiff = abs(numel(ref_autocorr_array) - numel(current_autocorr_array)); 
                        if (numel(ref_autocorr_array) < numel(current_autocorr_array))
                                ref_autocorr_array_padded = [ref_autocorr_array zeros(1, sizeDiff)]; % Zero Padding: Essential for KL Divergence
                                % Conversion To Probability Distribution: Essential for KL Divergence
                                prob_ref_autocorr = ref_autocorr_array_padded ./ sum(ref_autocorr_array_padded);
                                prob_current_autocorr = current_autocorr_array ./ sum(current_autocorr_array); 
                        else
                                current_autocorr_array_padded = [current_autocorr_array zeros(1, sizeDiff)]; % Zero Padding: Essential for KL Divergence
                                % Conversion To Probability Distribution: Essential for KL Divergence
                                prob_ref_autocorr = ref_autocorr_array ./ sum(ref_autocorr_array);
                                prob_current_autocorr = current_autocorr_array_padded ./ sum(current_autocorr_array_padded);
                        end
                            
                            % Compute the KL divergence of Autocorrelation
                            kl_divergence_autocorrel = sum(prob_ref_autocorr .* log2(prob_ref_autocorr ./ prob_current_autocorr), 'all');
                            normalized_divergence_autocorrel = 10*(kl_divergence_autocorrel / (numel(prob_ref_autocorr) + numel(prob_current_autocorr)));
                            rounded_divergence_autocorrel = round(normalized_divergence_autocorrel, 2);  % Round to 2 decimal places
                            % fprintf("****Normalized Autocorrelation function KL Divergence Score from Subject ID: %s Scenario: %s to Subject ID: %s Scenario: %s is: %f. \n \n", ID, scenario, current_ID_full, current_scenario, rounded_divergence_autocorrel);

                       %% Fundamental Frequency Section

                       current_fundamentalFrequency = interp1(ff_PS_w_norm, current_f, 1);
                       dev_fundamentalFrequency = abs(current_fundamentalFrequency - ref_fundamentalFrequency);

                       %% Second Breathing Harmonic Section

                       start_pos = 2 * round(current_fundamentalFrequency * 10) / 10;
                       end_pos = start_pos + 0.1;
                       breathHarmonic_freqRange_idx = (current_f >= start_pos) & (current_f <= end_pos);
                       breathHarmonic_values = interp1(current_f, ff_PS_w_norm, current_f(breathHarmonic_freqRange_idx));
                       [second_breathHarmonic, index] = max(breathHarmonic_values);
                       second_breathHarmonicFrequency = current_f(breathHarmonic_freqRange_idx);
                       current_second_breathHarmonicFrequency = second_breathHarmonicFrequency(index);
                       dev_second_breathHarmonicFrequency = abs(current_second_breathHarmonicFrequency - ref_second_breathHarmonicFrequency);

                       %% First Inter-Harmonic Distance Section

                       current_first_interHarmonicDist = current_second_breathHarmonicFrequency - current_fundamentalFrequency;
                       dev_first_interHarmonicDist = abs(current_first_interHarmonicDist - ref_first_interHarmonicDist);

                        %% Human Identification Prediction
                        diff = 0;
                        simi = 0;
                        if ((rounded_divergence_ps == 0) || (rounded_divergence_ps >= 0.07 && rounded_divergence_ps <= 0.16)) 
                            if ((dev_fundamentalFrequency >= 0.00 && dev_fundamentalFrequency <= 0.08) || (dev_second_breathHarmonicFrequency >= 0.00 && dev_second_breathHarmonicFrequency <= 0.03) || (dev_second_breathHarmonicFrequency >= 0.12 && dev_second_breathHarmonicFrequency <= 0.19) || (dev_second_breathHarmonicFrequency >= 0.21 && dev_second_breathHarmonicFrequency <= 0.22) || (dev_first_interHarmonicDist >= 0.01 && dev_first_interHarmonicDist <= 0.13))
                                if ((rounded_divergence_autocorrel >= -0.12 && rounded_divergence_autocorrel <= -0.11) || (rounded_divergence_autocorrel >= -0.02 && rounded_divergence_autocorrel <= 0.00) || (rounded_divergence_autocorrel == 0.02) || (rounded_divergence_autocorrel >= 0.06 && rounded_divergence_autocorrel <= 0.09) || ((rounded_divergence_autocorrel >= 0.22 && rounded_divergence_autocorrel <= 0.27)) || (rounded_divergence_autocorrel >= 0.52 && rounded_divergence_autocorrel <= 0.59))
                                    fprintf("-*-*-*-Prediction Result For Subject ID: %s Scenario: %s And Subject ID: %s Scenario: %s is: Same Subjects. \n \n", ID, scenario, current_ID_full, current_scenario);
                                    simi = 1;
                                else
                                    fprintf("-*-*-*-Prediction Result For Subject ID: %s Scenario: %s And Subject ID: %s Scenario: %s is: Different Subjects. \n \n", ID, scenario, current_ID_full, current_scenario);
                                    diff = 1;
                                end
                            else
                                fprintf("-*-*-*-Prediction Result For Subject ID: %s Scenario: %s And Subject ID: %s Scenario: %s is: Different Subjects. \n \n", ID, scenario, current_ID_full, current_scenario);
                                diff = 1;
                            end
                        else
                            fprintf("-*-*-*-Prediction Result For Subject ID: %s Scenario: %s And Subject ID: %s Scenario: %s is: Different Subjects. \n \n", ID, scenario, current_ID_full, current_scenario);
                            diff = 1;
                        end
                     
                        %% Confsuion Matrix Generation
                        % Similar is assumed to be Negative Class
                        % Different is assumed to be Positive Class
                        % TN is Ground Truth and Prediction are "Similar"
                        % FP is Ground Truth is "Similar" but Prediction is "Different"
                        % TP is Ground Truth and Prediction is "Different"
                        % FN is Ground Truth is "Positive" but Prediction is "Similar"

                        if (strcmp(current_ID_full, ID) == 0) % If scenario is different
                            if (diff == 1)
                                TP = TP + 1; % Ground Truth and Prediction are both Positive (Different)
                            elseif (simi == 1)
                                FN = FN + 1; % Ground Truth is Positive (Different), but Prediction is Negative (Similar)
                            end
                        else % If scenario is the same
                            if (simi == 1)
                                TN = TN + 1; % Ground Truth and Prediction are both Negative (Similar)
                            elseif (diff == 1)
                                FP = FP + 1; % Ground Truth is Negative (Similar), but Prediction is Positive (Different)
                            end
                        end

                end
            end
        end
        
        
    end
end

%% Confsuion Matrix Parametric Values
confusion_matrix = [TP, FP; FN, TN];
display(confusion_matrix);
sensitivity = (TP / (TP + FN))*100;
specificity = (TN / (FP + TN))*100;
precision = (TP / (TP + FP))*100;
neg_pred_value = (TN / (TN + FN))*100;
false_pos_rate = (FP / (FP + TN))*100;
false_discov_rate = (FP / (FP + TP))*100;
false_neg_rate = (FN / (FN + TP))*100;
accuracy =  ((TP + TN) / (TP + TN + FP +FN))*100;
f1_score = ((2*TP) / (2*TP + FP + FN))*100;

fprintf("Sensitivity of the Novel Human Identification Algorithm is: %f \n", sensitivity);
fprintf("Specificity of the Novel Human Identification Algorithm is: %f \n", specificity);
fprintf("Precision of the Novel Human Identification Algorithm is: %f \n", precision);
fprintf("Negative Predictive Value of the Novel Human Identification Algorithm is: %f \n", neg_pred_value);
fprintf("False Positive Rate of the Novel Human Identification Algorithm is: %f \n", false_pos_rate);
fprintf("False Discovery Rate of the Novel Human Identification Algorithm is: %f \n", false_discov_rate);
fprintf("False Negative Rate of the Novel Human Identification Algorithm is: %f \n", false_neg_rate);
fprintf("Accuracy of the Novel Human Identification Algorithm is: %f \n", accuracy);
fprintf("F1 Score of the Novel Human Identification Algorithm is: %f \n", f1_score);

% Plot confusion matrix
figure;
subplot(1, 2, 1);
imagesc(confusion_matrix);
colormap('redbluecmap');
colorbar('off');
textStrings = num2str(confusion_matrix(:), '%.1f');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:2);
hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));
textColors = repmat(confusion_matrix(:) > midValue, 1, 3);
set(hStrings, {'Color'}, num2cell(textColors, 2));
set(gca, 'XTick', 1:2, 'YTick', 1:2, 'XTickLabel', {'Different/Positive', 'Similar/Negative'}, 'YTickLabel', {'Different/Positive', 'Similar/Negative'}, 'FontSize', 10);
xlabel('True Label');
ylabel('Predicted Label');
title('Confusion Matrix');

% Display metrics as subplot texts
subplot(1, 2, 2);
textStrings = {sprintf('Accuracy: %.2f', accuracy), ...
               sprintf('Precision: %.2f', precision), ...
               sprintf('Sensitivity: %.2f', sensitivity), ...
               sprintf('Specificity: %.2f', specificity), ...
               sprintf('Neg Pred Value: %.2f', neg_pred_value), ...
               sprintf('False Pos Rate: %.2f', false_pos_rate), ...
               sprintf('False Discov Rate: %.2f', false_discov_rate), ...
               sprintf('False Neg Rate: %.2f', false_neg_rate), ...
               sprintf('F1 Score: %.2f', f1_score)};
text('Units', 'normalized', 'Position', [0.5, 0.5], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'String', textStrings);

% Adjust figure size and position
set(gcf, 'Position', [100, 100, 1200, 500]);
