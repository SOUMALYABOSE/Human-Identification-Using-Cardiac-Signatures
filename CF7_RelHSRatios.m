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
  IDrange = 18; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

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

        %% Cardiac Fingerprint 7: Relative Heart Sound Ratios
        figure;
        plot(time_radar,radar_heartsound.*10^6,'k-')
        hold on;
        set(vline(radar_hs_locsR(1)/fs_radar,'b-'),'HandleVisibility','on'); % Turn the legend on for vline
        vline(radar_hs_locsR./fs_radar,'b-');
        title('Compare heartbeat detection')
        xlim([0 2]);
        ylim([-200 200]);
        ylabel('Rel. Distance(um)');
        xlabel('Time(s)')
        legend('Radar heart sound','Radar S1')

        % Locating S-S Peak position
        radar_mod_heartsound = radar_heartsound.*10^6;  
        index_0pt4 = find(time_radar == 0.4, 1); % Find the index of 0.4 sec
        index_1pt0 = find(time_radar == 1.0, 1); % Find the index of 1.0 sec
        index_1pt5 = find(time_radar == 1.5, 1); % Find the index of 1.5 sec
        index_2pt0 = find(time_radar == 2, 1); % Find the index of 2.0 sec
        disp(['Index of 0.4 sec:', num2str(index_0pt4)]);
        disp(['Index of 1.0 sec:', num2str(index_1pt0)]);
        disp(['Index of 1.5 sec:', num2str(index_1pt5)]);
        disp(['Index of 2.0 sec:', num2str(index_2pt0)]);

        % First S-peak Value and Time-stamp
        start_index = index_0pt4;
        end_index = index_1pt0;
        search_range = radar_mod_heartsound(start_index:end_index);
        First_S_peak_value = min(search_range);
        disp('First S-peak Value:');
        disp(First_S_peak_value);
        index_First_S_peak = find(radar_mod_heartsound == First_S_peak_value, 1); % Find the index of First S-peak Value
        disp(['Index of First S-peak Value:', num2str(index_First_S_peak)]); % Index of First S-peak Value
        First_S_peak_timestamp = time_radar(index_First_S_peak);
        disp(['First S-peak Time-stamp:', num2str(First_S_peak_timestamp)]);

        % Second S-peak Value and Time-stamp
        start_index = index_1pt5;
        end_index = index_2pt0;
        search_range = radar_mod_heartsound(start_index:end_index);
        Second_S_peak_value = min(search_range);
        disp('Second S-peak Value:');
        disp(Second_S_peak_value);
        index_Second_S_peak = find(radar_mod_heartsound == Second_S_peak_value, 1); % Find the index of First S-peak Value
        disp(['Index of Second S-peak Value:', num2str(index_Second_S_peak)]); % Index of First S-peak Value
        Second_S_peak_timestamp = time_radar(index_Second_S_peak);
        disp(['Second S-peak Time-stamp:', num2str(Second_S_peak_timestamp)]);

        % Define the range of indices for peak detection
        start_index1 = index_First_S_peak;
        end_index1 = index_Second_S_peak;
        
        % Define the range of the array to analyze
        range_to_analyze = radar_mod_heartsound(start_index1:end_index1);
        
        % Find the peaks within the specified range
        [peaks, peak_indices] = findpeaks(range_to_analyze);
        
        % Adjust the peak indices to match the original array indices
        peak_indices = peak_indices + start_index1 - 1;
        
        % Plot the radar_mod_heartsound array
        figure;
        plot(time_radar,radar_mod_heartsound,'k-');
        hold on;
        set(vline(radar_hs_locsR(1)/fs_radar,'b-'),'HandleVisibility','on'); % Turn the legend on for vline
        vline(radar_hs_locsR./fs_radar,'b-');
        hold on;        
        %plot(time_radar(peak_indices), peaks, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); % Plot peaks as downward-pointing solid red triangles
        
        % Initialize arrays to store positive and negative peaks
        positive_peaks = peaks(peaks > 0);
        negative_peaks = peaks(peaks < 0);
        positive_peak_indices = peak_indices(peaks > 0);
        negative_peak_indices = peak_indices(peaks < 0);
        
        % Plot positive peaks as downward-pointing red triangles
        plot(time_radar(positive_peak_indices), positive_peaks, 'rv', 'MarkerSize', 10);
        
        % Plot negative peaks as upward-pointing red triangles
        plot(time_radar(negative_peak_indices), negative_peaks, 'r^', 'MarkerSize', 10);
        
        hold off;        
        % Label the peaks
        for i = 1:numel(peaks)
            text(time_radar(peak_indices(i)), peaks(i), ['Peak ' num2str(i)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        end        
        title('Peaks in radar_mod_heartsound')
        %xlim([time_radar(start_index1), time_radar(end_index1)]);
        xlim([0 2]);
        ylim([min(radar_mod_heartsound(start_index1:end_index1)), max(radar_mod_heartsound(start_index1:end_index1))]);
        ylabel('Rel. Distance(um)');
        xlabel('Time(s)')
        legend('Radar heart sound','Radar S1')
        % Display the plot
        grid on;
        

        
    end
end

