%% Important Copyright Information with regards to this scirpt:
% Part(s) of the Script has been taken from 2020 Script of: Sven Schellenberger, Kilin Shi (Copyright (C) 2020  Sven Schellenberger, Kilin Shi)
% Rest has been added under the CME-Master Thesis titled "Human Identification With Gender Knowledge Analysis From Cardiac Signatures and Health-related Information Extraction From Radar Based Biomarkers" 
% Conducted @ Home Automation Lab, Institute FAPS, FAU Erlangen-Nuremberg by Mr. Soumalya Bose
% Advisors: Mr. Jochen Bauer (Institute FAPS, FAU), Dr. -med. Tobias Steigleder (PallMeT, Univ. Hospital Erlangen, FAU), Professor Dr.-Ing. Jörg Franke (Institute FAPS, FAU), Professor Dr.-Ing. Georg Fischer (Institute for Electronics Engineering, FAU)
% Copyright (C) 2023-24  Mr. Soumalya Bose
%% Clearing of Workspace
clear;
clc;

%% Init
addpath(genpath('utils'))


%% Preparation of Datasets

%%%% Or, Possibility 2:  Choose Subject-ID(s)  (Uncomment/Comment When Needed)
% GDN00XX - simply choose numbers or ranges from 01 to 30
% IDrange = [2 3 4 5 6 7 8 9 10 12 13 15]; % Sample (Uncomment/Comment When Needed)
  IDrange = 30; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
% scenarios = {'Resting', 'Valsalva', 'Apnea'};
  scenarios = {'Resting'}; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
 

%%%% Set path to datasets 
% Datasets can be found on figshare
path = 'datasets'; % In this case "datasets" folder is in this scripts folder 
scrsz = get(groot,'ScreenSize'); % For plotting

%% Extracting Info from Datasets with further Processing and Ploting of relevant Data 

for indx = 1:length(IDrange)
    %% Iterate over all subject IDs
    ID = sprintf('GDN%04d',IDrange(indx));
    fprintf('----------- Loading %s ------------\n', ID);

    for sz = 1:length(scenarios)
        %% Iterate over all existing subject scenarios for the selected subject ID
        scenario = scenarios{sz};
        fprintf('---- Scenario %s\n', scenario);

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
            ref_PS_w_norm = 0; % Customized by Soumalya Bose to avoid Technical Glitches
            continue
        end

        %% Prepare data
        output.(ID).(scenario) = struct;
        [~,~,~,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},0); % Ellipse fitting, compensation and distance reconstruction

        % Radar
        [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);

        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar; 

        %% Pre-processing of Radar Data/Rel. Distance/Chest Displacement Signal 
        radar_dist_scaled = radar_dist*1000; % Scaling of Raw Radar Displacement Signal
        % Kaiser Window with β = 5.1 (0-4: Wide ML-L. Suppre SL-High Freq Resol. 4-10: Reverse+Less Spectral Leak/Artefact)
        w = kaiser(length(radar_dist_scaled), 5.1); 
        radar_dist_windowed = radar_dist_scaled.*w;
        Yw = fft(radar_dist_windowed);
        L = length(radar_dist_windowed);

        %% Fequency Domain Transformation of Processed Radar Data/Rel. Distance/Chest Displacement Signal 
        PS2w = abs(Yw/L); % Double sampling plot
        PS1w = PS2w(1:floor(L/2)+1); % Single sampling plot. Used Floor func. as L is from length func. which returns double and not integer.
        PS1w(2:end-1) = 2*PS1w(2:end-1);
        f = fs_radar*(0:(L/2))/L;  % Frequency bins

        %% Normalization of Magnitude Spectrum from 0.1 Hz to 1.0 Hz
        PS1_vals = interp1(f,PS1w,0.1:0.00004:1);
        max_PS1w = max(PS1_vals);
        PS_w_norm = PS1w/max_PS1w;
        ref_PS_w_norm = interp1(f,PS_w_norm,0.1:0.00004:0.8);

        %% Cardiac Fingerprint 5: Kullback-Leibler (KL) Divergence Scores (Applicable for Human Identification Analysis)
        % Divergence Scores are computed with Reference Scenario: Resting only in consideration. This Scenarios has the least Radar Data noise
        % Measures difference between two Probability distributions in terms of their information content
        % Larger values indicate greater divergence 
        % Smaller values suggest a higher degree of similarity

            if (scenario == "Resting")
                current_IDrange = IDrange; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
                current_scenarios = {'Valsalva', 'Apnea'};

                for current_indx = 1:length(current_IDrange)
                    current_ID = current_IDrange(current_indx);
                    current_ID_full = sprintf('GDN%04d',current_ID);
                    for current_sz = 1:length(current_scenarios)
                        current_scenario = current_scenarios{current_sz};
                        current_PS_w_norm = Ps_w_norm_calc(current_ID, current_scenario);

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
                        % Compute the KL divergence
                        kl_divergence = sum(prob_ref_PS_w_norm .* log2(prob_ref_PS_w_norm ./ prob_current_PS_w_norm), 'all');
                        normalized_divergence = 10000*(kl_divergence / (numel(prob_ref_PS_w_norm) + numel(prob_current_PS_w_norm)));
                        rounded_divergence = round(normalized_divergence, 2);  % Round to 2 decimal places
                        fprintf("****Normalized KL Divergence Score from Subject ID: %s Scenario: %s to Subject ID: %s Scenario: %s is: %f. \n", ID, scenario, current_ID_full, current_scenario, rounded_divergence);
                    end
                end
            end

    end
 end