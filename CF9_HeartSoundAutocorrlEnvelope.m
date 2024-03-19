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
  IDrange = 12; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

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

        %% Cardiac Fingerprint 7: Heart Sound Time Autocorrelation Function Envelope
        radar_mod_heartsound = radar_heartsound.*10^6;

        % Calculate autocorrelation
        [autocorr_values_plot, lag_indices_plot] = xcorr(radar_mod_heartsound, 'biased');
        
        % Plot autocorrelation
        figure;
        plot(lag_indices_plot, autocorr_values_plot);
        title('Autocorrelation of Radar Heart Sound');
        xlim([0 2000]);
        xlabel('Lag');
        ylabel('Autocorrelation');
        hold on;
        
        % % Calculate envelope of autocorrelation
        % [upper_env, lower_env] = envelope(autocorr_values_plot, 100, 'peak');
        % hold on;
        % plot(lag_indices_plot, upper_env, 'r', 'LineWidth', 2);
        % plot(lag_indices_plot, lower_env, 'r', 'LineWidth', 2);
        % legend('Autocorrelation', 'Envelope');
        % hold off;
        % 
        % % Fit plot to screen
        % axis tight;
        % 
        % % Calculate area under upper envelope using trapezoidal rule
        % area_upper = trapz(lag_indices_plot, upper_env);
        % 
        % % Calculate area under lower envelope using trapezoidal rule
        % area_lower = trapz(lag_indices_plot, lower_env);
        % 
        % Total area enclosed by envelope
        % total_area = abs(area_upper - area_lower);
        % disp(['Total area enclosed by envelope: ', num2str(total_area)]);
        
        
    end
end

