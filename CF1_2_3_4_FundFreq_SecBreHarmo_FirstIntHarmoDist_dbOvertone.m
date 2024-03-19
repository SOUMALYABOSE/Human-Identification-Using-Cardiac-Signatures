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
  IDrange = 28; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
% scenarios = {'Resting', 'Valsalva', 'Apnea', 'TiltUp', 'TiltDown'}; % (Uncomment/Comment When Needed)
  scenarios = {'Valsalva'}; % 1-2-3 Testing Scenario. (Uncomment/Comment When Needed)
 

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

        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar; 

        %% Pre-processing of Radar Data/Rel. Distance/Chest Displacement Signal 
        radar_dist_scaled = radar_dist*1000; % Scaling of Raw Radar Displacement Signal
        % Kaiser Window with β = 5.1 (0-4:Wide ML-L.Suppre SL-High Freq Resol. 4-10:Reverse+Less Spectral Leak/Artefact
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
         
        %% Plotting of Normalized Magnitude Spectrum 
        figure;
        plot(f, PS_w_norm);
        box_color = 'red';
        box_opacity = 0.3;
        annotation('textbox', [0.8, 0.8, 0.1, 0.1], 'String', scenario, 'Color', 'black', 'EdgeColor', box_color, 'BackgroundColor', box_color, 'FaceAlpha', box_opacity, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        xlim([0.1 1]); % Breathing: 0.2 Hz means 12 (60*0.2=12) . Heartbeat: 1 Hz
        ylim([0 1]);
        title("Single Sided Mag. Spectrum With Kaiser Windowing (β = 5.1) for Subject ID:", ID);
        xlabel("Frequency (Hz)");
        ylabel("Normalized Magnitude Spectrum 0.1 Hz-1.0 Hz");

        %% Cardiac Fingerprint 1: Fundamental Frequency (Applicable for Gender Analysis and Human Identification Analysis)
        fundamentalFrequency = interp1(PS_w_norm, f, 1);
        fprintf("****Fundamental Frequency Position is: %f Hz. \n", fundamentalFrequency);

        %% Cardiac Fingerprint 2: Second Breathing Harmonic (Applicable for Gender Analysis and Human Identification Analysis)
        start_pos = 2 * round(fundamentalFrequency * 10) / 10;
        end_pos = start_pos + 0.1;
        breathHarmonic_freqRange_idx = (f >= start_pos) & (f <= end_pos);
        breathHarmonic_values = interp1(f, PS_w_norm, f(breathHarmonic_freqRange_idx));
        [second_breathHarmonic, index] = max(breathHarmonic_values);
        second_breathHarmonicFrequency = f(breathHarmonic_freqRange_idx);
        second_breathHarmonicFrequency = second_breathHarmonicFrequency(index);
        fprintf("****2nd Breathing Harmonic Position is: %f Hz. \n", second_breathHarmonicFrequency);

        %% Cardiac Fingerprint 3: First Inter-Harmonic Distance (Applicable for Gender Analysis and Human Identification Analysis)
        first_interHarmonicDist = second_breathHarmonicFrequency - fundamentalFrequency;
        fprintf("****1st Inter-Harmonic Distance is: %f Hz. \n", first_interHarmonicDist);
        
        %% Cardiac Fingerprint 4: dB Overtone (Applicable for Gender Analysis)
        db_overtone = 20*log10(1 - second_breathHarmonic);
        fprintf("****dB Overtone is: %f Hz. \n", db_overtone); % Likely a Negative value as 2nd Breathing Harmonic is always less than reference value i.e. Fundamental Frequency
      % fprintf("****dB Overtone for Subject ID: %s in Scenario: %s is: %f Hz. \n",ID,scenario,db_overtone); % Sample

    end
end
