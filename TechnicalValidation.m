%% Authors: Sven Schellenberger, Kilin Shi
% Copyright (C) 2020  Sven Schellenberger, Kilin Shi

%%%%%% IMPORTANT %%%%%%
% Copy the data from figshare (link to be generated) to a subfolder "datasets" inside the folder of
% this script or change datadir to according dataset location

% Set the subject-ID(s) and scneraio(s) you want to view in the
% 'Options'-section and start the script
%%%%%%%%%%%%%%%%%%%%%%%

%% Validation
% This script builds a database of interbeat intervals (IBIs) evaluated
% from radar heart sound and reference ecg (first section)
% and evaluates the data of all subjects in the second sections
% The plots seen in the publication are generated.
% If database already exists, creation will be skipped.

%% Init
addpath(genpath('utils'))

%% Options
%%%% Choose subject-ID(s) 
IDrange = 1:30; % GDN00XX, simply choose numbers or ranges from 1 to 30
% IDrange = 14; % GDN00XX

%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}
scenarios = {'Resting'}; % Choose one or more scenarios for evaluation 

%%%% Set path to datasets 
% Datasets can be found on figshare
path = 'datasets'; % In this case "datasets" folder is in this scripts folder 

%% Build database 
% Database construction for an easier evaluation of data

% Manually selected ECG channel for each subject (Position in ECG_CHANNEL vector corresponds to ID)
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
output = struct;

if ~exist('output_validation.mat','file')
for indx = 1:length(IDrange)
    %% Iterate all subject IDs
    ID = sprintf('GDN%04d',IDrange(indx));
    fprintf('----------- Loading %s ------------\n', ID);
    output.(ID) = struct;
    
    for sz = 1:length(scenarios)
        %% Iterate all existing subject scenarios
        scenario = scenarios{sz};
        fprintf('--- %s\n', scenario);
        
        % Search file
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
        
        % Reconstruct raw radar data to fit on unit circle
        [~,~,~,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},0);
        % Usage of elreko
        % [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1}(timestamp of dataset),0(Plot flag -> 1: plots on, 0: plots off));
        
        % In scenarios TiltUp and TiltDown the tilt movement needs to be
        % removed for evaluation
        sc_sec = getInterventions(scenario,tfm_intervention,fs_intervention);
        if strcmp(scenario,'TiltUp') || strcmp(scenario,'TiltDown') 
            % cut fs == 2000
            radar_dist = radar_dist(sc_sec(1):sc_sec(2));
            tfm_ecg1 = tfm_ecg1(sc_sec(1):sc_sec(2));
            tfm_ecg2 = tfm_ecg2(sc_sec(1):sc_sec(2));
            % cut fs == 1000
            if abs(floor(sc_sec(2)/2) - length(tfm_icg)) < 10
                tfm_icg = tfm_icg(ceil(sc_sec(1)/2):length(tfm_icg));
            else
                tfm_icg = tfm_icg(ceil(sc_sec(1)/2):floor(sc_sec(2)/2));
            end
            % cut fs == 200
            if abs(floor(sc_sec(2)/10) - length(tfm_bp)) < 10
                tfm_bp = tfm_bp(ceil(sc_sec(1)/10):length(tfm_bp));
            else
                tfm_bp = tfm_bp(ceil(sc_sec(1)/10):floor(sc_sec(2)/10));
            end
            % cut fs == 100
            if abs(floor(sc_sec(2)/20) - length(tfm_z0)) < 10
                tfm_z0 = tfm_z0(ceil(sc_sec(1)/20):length(tfm_z0));
            else
                tfm_z0 = tfm_z0(ceil(sc_sec(1)/20):floor(sc_sec(2)/20));
            end
        end
        
        [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);

        % Process ECG signals
        if ECG_CHANNEL(IDrange(indx)) == 1
            tfm_ecg = fillmissing(tfm_ecg1,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        else
            tfm_ecg = fillmissing(tfm_ecg2,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        end
        tfm_ecg = filtButter(tfm_ecg,fs_ecg,4,[1 20],'bandpass');
        
        % Resample radar respiration to match tfm respiration
        radar_respiration_re = resample(radar_respiration,fs_z0,fs_radar);
        radar_dist_re = resample(radar_dist,fs_z0,fs_radar);
        if length(radar_respiration_re) > length(tfm_respiration)
            radar_respiration_re = radar_respiration_re(1:length(tfm_respiration));
            radar_dist_re = radar_dist_re(1:length(tfm_respiration));
        elseif length(radar_respiration_re) < length(tfm_respiration)
            tfm_respiration = tfm_respiration(1:length(radar_respiration_re));
        end
        
        
        %% Evaluate data
        % Heartbeat detection
        % Determine heart sound in radar signal and extract S1 location 
        radar_hs_states = getHsmmStates(radar_heartsound, fs_radar);
        radar_hs_locsR = statesToLocs(radar_hs_states, 1);
        % Detect R-peaks in ECG signal
        [~,tfm_ecg_locsR] = twaveend(tfm_ecg(1:length(radar_hs_states)), fs_ecg,32*(fs_ecg/200),'p');
        % Save subject data 
        [output.(ID).(scenario).tfm_ecg_ibi, output.(ID).(scenario).radar_hs_ibi] = getIbisNoInterp(tfm_ecg_locsR, radar_hs_locsR, 5, 6, fs_radar);
        fprintf('---- done\n');
    end
end
    save('output_validation.mat','output'); 
else
    load('output_validation.mat');
end

%% Output data
% choose other scenarios for evaluation
% scenarios = {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'};
% scenarios = {'Resting'}; % 

persons = fieldnames(output);
for sz = 1:length(scenarios)
    %% Gather precalculated data
    scenario = scenarios{sz};
    tfm_ecg_ibi_overall = [];
    radar_hs_ibi_overall = [];
    rmse_per_person = [];
    for ct_persons = 1:numel(persons)
        person = char(persons(ct_persons));
        if isfield(output.(person),scenario)
            tfm_ecg_ibi_overall = [tfm_ecg_ibi_overall; output.(person).(scenario).tfm_ecg_ibi];
            radar_hs_ibi_overall = [radar_hs_ibi_overall; output.(person).(scenario).radar_hs_ibi];
            diff_radar_ecg_person = abs(output.(person).(scenario).tfm_ecg_ibi - output.(person).(scenario).radar_hs_ibi);
            rmse_per_person = [rmse_per_person; sqrt(mean(diff_radar_ecg_person.^2))*1000];
        end
    end
    
    %% Evaluate gathered data
    % Data as seen in publication
    
    corr_radar_ecg = corr(tfm_ecg_ibi_overall,radar_hs_ibi_overall);
    diff_radar_ecg = abs(tfm_ecg_ibi_overall - radar_hs_ibi_overall);
    rmse_radar_ecg = sqrt(mean(diff_radar_ecg.^2))*1000;
    
    % Bland-Altman
    [means,diffs,meanDiff,CR,~] = BlandAltman(radar_hs_ibi_overall',tfm_ecg_ibi_overall',2,0);     

    % Histogram
    figure;hnd = histogram(diff_radar_ecg);

    % RMSE
    figure;bar(rmse_per_person);
    hold on;
    hline(median(rmse_per_person));

    % Scatterplot
    selgreen = diff_radar_ecg < 0.01;
    selblue = diff_radar_ecg >= 0.01 & diff_radar_ecg < 0.075;
    selred = diff_radar_ecg > 0.075;

    figure;
    scatter(tfm_ecg_ibi_overall(selgreen)',radar_hs_ibi_overall(selgreen)',2,'filled','g');
    hold on;
    scatter(tfm_ecg_ibi_overall(selblue)',radar_hs_ibi_overall(selblue)',2,'filled','b');
    scatter(tfm_ecg_ibi_overall(selred)',radar_hs_ibi_overall(selred)',2,'filled','r');
    diagonal_line = 0.6:0.1:1.6;
    plot(diagonal_line,diagonal_line);
    
    str_gr = sprintf('IBIdiff < 0.01s : %.1f%%', (sum(selgreen)/length(diff_radar_ecg))*100);
    str_bl = sprintf('0.01s < IBIdiff < 0.075s : %.1f%%', (sum(selblue)/length(diff_radar_ecg)*100));
    str_re = sprintf('0.075s IBIdiff < 0.01s : %.1f%%', (sum(selred)/length(diff_radar_ecg)*100));
    legend(str_gr,str_bl,str_re,'Diagonal','Location','southeast')
end
