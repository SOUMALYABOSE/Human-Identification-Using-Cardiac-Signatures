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
% IDrange = [1 11 14 17 19 21 22 29 30]; (Uncomment/Comment When Needed)
  IDrange = 28; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
% scenarios = {'Resting', 'Valsalva', 'Apnea', 'TiltUp', 'TiltDown'}; % (Uncomment/Comment When Needed)
  scenarios = {'Resting'}; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
 

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
    fprintf('----------- Loading %s ------------\n', ID);

    for sz = 1:length(scenarios)
        %% Iterate all existing subject scenarios
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
            ref_autocorr_array = 0; % Customized by Soumalya Bose to avoid Technical Glitches
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

        %% Cardiac Fingerprint 7: Heart Sound Time Axis Analysis
        radar_mod_heartsound = radar_heartsound.*10^6;
        % figure;
        % plot(time_radar,radar_mod_heartsound,'k-')
        % hold on;
        % set(vline(radar_hs_locsR(1)/fs_radar,'b-'),'HandleVisibility','on'); % Turn the legend on for vline
        % vline(radar_hs_locsR./fs_radar,'b-');
        % title('Compare heartbeat detection')
        % xlim([0 2]);
        % ylim([-200 200]);
        % ylabel('Rel. Distance(um)');
        % xlabel('Time(s)')
        % legend('Radar heart sound','Radar S1')

        % Plot autocorrelation function
        % figure;
        % [autocorr_values_plot, lag_indices_plot] = xcorr(radar_mod_heartsound, 'biased');
        % plot(lag_indices_plot, autocorr_values_plot);
        % title('Autocorrelation of Radar Heart Sound');
        % xlabel('Lag');
        % ylabel('Autocorrelation');

        % Calculate autocorrelation
        lags = -length(radar_mod_heartsound) + 1 : length(radar_mod_heartsound) - 1;
        start_index = 1;  % Specify the starting index for autocorrelation calculation
        end_index = min(start_index + 10000, length(radar_mod_heartsound));  % Specify the ending index, limit to a reasonable size
        [autocorr_values, lag_indices] = xcorr(radar_mod_heartsound(start_index:end_index), 'biased');
        ref_autocorr_array = autocorr_values;

        % Calculating KL Divergence between Different Autocorrelation Functions cases
        if (scenario == "Resting")
            current_IDrange = [5 8 12 28]; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
            current_scenarios = {'Resting', 'Valsalva', 'Apnea'};

            for current_indx = 1:length(current_IDrange)
                current_ID = current_IDrange(current_indx);
                current_ID_full = sprintf('GDN%04d',current_ID);
                for current_sz = 1:length(current_scenarios)
                    current_scenario = current_scenarios{current_sz};
                    current_autocorr_array = Autocorrl_calc(current_ID, current_scenario);

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
                        % Compute the KL divergence
                        kl_divergence = sum(prob_ref_autocorr .* log2(prob_ref_autocorr ./ prob_current_autocorr), 'all');
                        normalized_divergence = 10*(kl_divergence / (numel(prob_ref_autocorr) + numel(prob_current_autocorr)));
                        rounded_divergence = round(normalized_divergence, 2);  % Round to 2 decimal places
                        fprintf("****Normalized Autocorrelation function KL Divergence Score from Subject ID: %s Scenario: %s to Subject ID: %s Scenario: %s is: %f. \n", ID, scenario, current_ID_full, current_scenario, rounded_divergence);

                end
            end
        end
        
        
    end
end

