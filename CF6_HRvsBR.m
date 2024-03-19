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

% GDN00XX - simply choose numbers or ranges from 01 to 30
% IDrange = [1 11 14 17 19 21 22 29 30]; (Uncomment/Comment When Needed)
  IDrange = 3; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
% scenarios = {'Resting', 'Valsalva', 'Apnea', 'TiltUp', 'TiltDown'}; % (Uncomment/Comment When Needed)
  scenarios = {'Resting', 'Valsalva', 'Apnea'}; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
 

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
            continue
        end

        %% Prepare data
        output.(ID).(scenario) = struct;
        [~,~,~,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},0); % Ellipse fitting, compensation and distance reconstruction

        % Radar
        [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);

        %% TFM
        % Manually selected ECG channel for each subject (Position in ECG_CHANNEL vector corresponds to ID)
        ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];

        if ECG_CHANNEL(IDrange(indx)) == 1
            tfm_ecg = fillmissing(tfm_ecg1,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        else
            tfm_ecg = fillmissing(tfm_ecg2,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        end
        tfm_ecg = filtButter(tfm_ecg,fs_ecg,4,[1 20],'bandpass');
        
        %% Resample radar respiration to match tfm respiration
        radar_respiration_re = resample(radar_respiration,fs_z0,fs_radar);
        radar_dist_re = resample(radar_dist,fs_z0,fs_radar);
        if length(radar_respiration_re) > length(tfm_respiration)
            radar_respiration_re = radar_respiration_re(1:length(tfm_respiration));
            radar_dist_re = radar_dist_re(1:length(tfm_respiration));
        elseif length(radar_respiration_re) < length(tfm_respiration)
            tfm_respiration = tfm_respiration(1:length(radar_respiration_re));
        end
        sc_sec = getInterventions(scenario,tfm_intervention,fs_intervention);

        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar; 
        time_ecg = 1/fs_ecg:1/fs_ecg:length(tfm_ecg)/fs_ecg;
        time_icg = 1/fs_icg:1/fs_icg:length(tfm_icg)/fs_icg;
        time_bp = 1/fs_bp:1/fs_bp:length(tfm_bp)/fs_bp;
        time_z0 = 1/fs_z0:1/fs_z0:length(tfm_z0)/fs_z0;
        time_respiration = 1/fs_z0:1/fs_z0:length(radar_respiration_re)/fs_z0;

        %% Cardiac Fingerprint 5: Average Heart Rate (Calculated from Heart Sounds)
        radar_hs_states = getHsmmStates(radar_heartsound, fs_radar);
        radar_hs_locsR = statesToLocs(radar_hs_states, 1);

        radar_heartrate = 60./(diff(radar_hs_locsR)./fs_radar);
        radar_heartrate_mean = mean(radar_heartrate, 'all'); % 40-150 BPM Expected for Subject IDs (Healthy)
        fprintf("****Average Heart Rate (BPM) for Subject ID: %s in Scenario: %s is: %f. \n",ID,scenario,radar_heartrate_mean);

        figure;
        subplot(2,1,1);
        plot(radar_hs_locsR(1:end-1)./fs_radar,60./(diff(radar_hs_locsR)./fs_radar));
        title('Heart rate');
        xlabel('Time(s)');
        ylabel('Heart rate(BPM)');
        hold on;

        % Add text annotation for average heart rate
        text(0.5, 0.95, sprintf('Average Heart Rate: %.2f BPM', radar_heartrate_mean), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'normalized');

        %% Cardiac Fingerprint 6: Average Respiration Rate
        radar_mod_respiration = radar_respiration_re.*1000;
   
        subplot(2,1,2);
        plot(time_respiration, radar_respiration_re.*1000, '-k');
        xlim([0 60]);
        title('Compare respiration');
        xlabel('Time(s)');
        ylabel('Rel. Distance(mm)');
        hold on;

        % Calculate time_range_indices
        time_range_indices = find(time_respiration >= 0 & time_respiration <= 60);
        
        % Plot radar_mod_respiration
        plot(time_respiration(time_range_indices), radar_mod_respiration(time_range_indices), '-b');
        
        % Find peaks within the specified time range
        [peaks, peak_indices] = findpeaks(radar_mod_respiration(time_range_indices));

        fprintf("****Average Breathing Rate per Min. for Subject ID: %s in Scenario: %s is: %f. \n",ID,scenario,numel(peak_indices));
        
        % Plot the peaks on the existing plot
        plot(time_respiration(time_range_indices(peak_indices)), peaks, 'ro', 'MarkerFaceColor', 'r');

        % Label the peaks with numbers
        for i = 1:numel(peak_indices)
            text(time_respiration(time_range_indices(peak_indices(i))), peaks(i), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end

        % Annotate the plot
        legend('','Radar Respiration', 'Peaks');
        hold off;
        

    end
end
