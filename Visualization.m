%% Important Copyright Information with regards to this scirpt:
% Part(s) of the Script has been taken from 2020 Script of: Sven Schellenberger, Kilin Shi (Copyright (C) 2020  Sven Schellenberger, Kilin Shi)
% Rest has been added under the CME-Master Thesis titled "Human Identification With Gender Knowledge Analysis From Cardiac Signatures and Health-related Information Extraction From Radar Based Biomarkers" 
% Conducted @ Home Automation Lab, Institute FAPS, FAU Erlangen-Nuremberg by Mr. Soumalya Bose
% Advisors: Mr. Jochen Bauer (Institute FAPS, FAU), Dr. -med. Tobias Steigleder (PallMeT, Univ. Hospital Erlangen, FAU), Professor Dr.-Ing. Jörg Franke (Institute FAPS, FAU), Professor Dr.-Ing. Georg Fischer (Institute for Electronics Engineering, FAU)
% Copyright (C) 2023-24  Mr. Soumalya Bose
%% Clearing of Workspace
clear;
clc;

%% Declaring all arrays for Human Identification Related Visualization

% Subject IDs
subjects = {'GDN0001', 'GDN0002', 'GDN0003', 'GDN0004', 'GDN0006', 'GDN0007', 'GDN0009', 'GDN0010', 'GDN0011', 'GDN0013', 'GDN0014', 'GDN0015', 'GDN0018', 'GDN0020', 'GDN0021', 'GDN0022', 'GDN0024', 'GDN0025', 'GDN0026', 'GDN0027', 'GDN0029'};
subjects_full = {'GDN0001', 'GDN0001', 'GDN0001', 'GDN0002', 'GDN0002','GDN0002', 'GDN0003', 'GDN0003','GDN0003', 'GDN0004', 'GDN0004', 'GDN0004', 'GDN0006', 'GDN0006', 'GDN0006', 'GDN0007', 'GDN0007', 'GDN0007', 'GDN0009', 'GDN0009', 'GDN0009', 'GDN0010', 'GDN0010', 'GDN0010', 'GDN0011', 'GDN0011', 'GDN0011', 'GDN0013', 'GDN0013', 'GDN0013', 'GDN0014', 'GDN0014', 'GDN0014', 'GDN0015', 'GDN0015', 'GDN0015', 'GDN0018', 'GDN0018','GDN0018', 'GDN0020', 'GDN0020','GDN0020', 'GDN0021', 'GDN0021','GDN0021', 'GDN0022', 'GDN0022','GDN0022', 'GDN0024', 'GDN0024','GDN0024', 'GDN0025', 'GDN0025','GDN0025', 'GDN0026', 'GDN0026','GDN0026', 'GDN0027', 'GDN0027','GDN0027', 'GDN0029', 'GDN0029', 'GDN0029'};


% Fundamental frequencies Info
fundamental_frequency = [0 0 0 0 0.06 0 0.01 0.02 0.13 0.07 0.04 0 0.07 0.02 0.02 0.08 0 0.06 0 0.11 0.01];

% Second Breathing Harmonic Final Relative Difference
second_breathing_harmonic = [0 0 0 0.16 0 0.01 0.01 0.01 0.12 0.19 0.21 0 0.21 0.03 0.19 0.25 0 0.07 0 0.22 0.15];

% First Inter-Harmonic Distance Final Relative Difference
first_interHarmonic_distance = [0 0 0 0.13 0.06 0.02 0.03 0.02 0.06 0.10 0 0.12 0.05 0.13 0.12 0.01 0 0.05 0.22 0.13 0.07];

% Power Spectrum
ps = [0.0 0.09 0 0.0 0.07 0 0.0	0.11 0 0.0 0.14	0.16 0.0 0.14 0.16 0.0 0.09 0.13 0.0 0.11 0.13 0.0 0.13 0.25 0.0 0.16 0.16 0.0 0.10 0.10 0.0 0.12 0.15 0.0 0 0 0.0 0.10 0.08 0.0 0.10 0.10 0.0 0.17 0.35 0.0 0.14 0.26 0.0 0 0 0.0 0.09 0.10 0.0 0 0 0.0 0.15 0.08 0.0 0.18 0.09];
 
%% Visualization Fundamental Frequency for Human Identification
threshold_min_ff = 0.0;
threshold_max_ff = 0.08;

% Create scatter plot
figure;
scatter(1:length(subjects), fundamental_frequency, 'filled');

% Set x-axis labels
xticks(1:length(subjects));
xticklabels(subjects);
xlabel('Subject ID');

% Set y-axis label
ylabel('Fundamental Frequency (Hz)');

% Add horizontal dotted threshold lines
hold on;
line([0, length(subjects)+1], [threshold_min_ff, threshold_min_ff], 'LineStyle', '--', 'Color', 'k');
line([0, length(subjects)+1], [threshold_max_ff, threshold_max_ff], 'LineStyle', '--', 'Color', 'k');
hold off;

% Adjust y-axis limits if necessary
ylim([min(fundamental_frequency)-0.1, max(fundamental_frequency)+0.1]);

% Add legend
legend('Fundamental Frequency', 'Lower Threshold', 'Upper Threshold');

% Add grid
grid on;

% Title
title('Fundamental Frequency Relative Difference Scatter Plot for Human Identification');

%% Visualization Second Breathing Harmonic for Human Identification
threshold1_sbh = 0.0;
threshold2_sbh = 0.03;
threshold3_sbh = 0.12;
threshold4_sbh = 0.19;
threshold5_sbh = 0.21;

% Create scatter plot
figure;
scatter(1:length(subjects), second_breathing_harmonic, 'filled');
hold on; % Allow multiple plots on the same axes

% Set x-axis labels
xticks(1:length(subjects));
xticklabels(subjects);
xlabel('Subject ID');

% Set y-axis label
ylabel('Second Breathing Harmonic (Hz)');

% Adjust y-axis limits if necessary
ylim([min(second_breathing_harmonic)-0.1, max(second_breathing_harmonic)+0.1]);

% Add grid
grid on;

% Plot threshold lines
line([1, length(subjects)], [threshold1_sbh, threshold1_sbh], 'Color', 'r', 'LineStyle', '--');
line([1, length(subjects)], [threshold2_sbh, threshold2_sbh], 'Color', 'r', 'LineStyle', '--');
line([1, length(subjects)], [threshold3_sbh, threshold3_sbh], 'Color', 'g', 'LineStyle', '--');
line([1, length(subjects)], [threshold4_sbh, threshold4_sbh], 'Color', 'g', 'LineStyle', '--');
line([1, length(subjects)], [threshold5_sbh, threshold5_sbh], 'Color', 'b', 'LineStyle', '--');

% Add legend for threshold lines with custom labels
legend('Second Breathing Harmonic Relative Difference', 'Threshold 1', 'Threshold 2', 'Threshold 3', 'Threshold 4', 'Threshold 5', 'Location', 'NorthEast');

% Title
title('Second Breathing Harmonic Relative Difference Scatter Plot for Human Identification');

%% Visualization First Inter-Harmonic Distance for Human Identification
threshold1_fih = 0.0;
threshold2_fih = 0.13;

% Create scatter plot
figure;
scatter(1:length(subjects), first_interHarmonic_distance, 'filled');
hold on; % Allow multiple plots on the same axes

% Set x-axis labels
xticks(1:length(subjects));
xticklabels(subjects);
xlabel('Subject ID');

% Set y-axis label
ylabel('First Inter-Harmonic Distance (Hz)');

% Adjust y-axis limits if necessary
ylim([min(first_interHarmonic_distance)-0.1, max(first_interHarmonic_distance)+0.1]);

% Add grid
grid on;

% Plot threshold lines
line([1, length(subjects)], [threshold1_fih, threshold1_fih], 'Color', 'r', 'LineStyle', '--');
line([1, length(subjects)], [threshold2_fih, threshold2_fih], 'Color', 'r', 'LineStyle', '--');

% Add legend for threshold lines with custom labels
legend('First Inter-Harmonic Distance Relative Difference', 'Threshold 1', 'Threshold 2', 'Location', 'NorthEast');

% Title
title('First Inter-Harmonic Distance Relative Difference Scatter Plot for Human Identification');

%% Visualization Power Spectrum KL Divergence for Human Identification
threshold1_ps = 0.0;
threshold2_ps = 0.07;
threshold3_ps = 0.16;

% Create scatter plot
figure;
scatter(1:length(subjects_full), ps, 'filled');
hold on; % Allow multiple plots on the same axes

% Set x-axis labels
xticks(1:length(subjects_full));
xticklabels(subjects_full);
xlabel('Subject ID');

% Set y-axis label
ylabel('Power Spectrum KL Divergence Score');

% Adjust y-axis limits if necessary
ylim([min(ps)-0.1, max(ps)+0.1]);

% Add grid
grid on;

% Plot threshold lines
line([1, length(subjects_full)], [threshold1_ps, threshold1_ps], 'Color', 'r', 'LineStyle', '--');
line([1, length(subjects_full)], [threshold2_ps, threshold2_ps], 'Color', 'b', 'LineStyle', '--');
line([1, length(subjects_full)], [threshold3_ps, threshold3_ps], 'Color', 'b', 'LineStyle', '--');

% Add legend for threshold lines with custom labels
legend('Power Spectrum KL Divergence Score', 'Threshold 1', 'Threshold 2', 'Threshold 3', 'Location', 'NorthEast');

% Title
title('Power Spectrum KL Divergence Score Scatter Plot for Human Identification');

%% Declaring all arrays for Gender Classification Analysis Related Visualization
subjects_gender = {'GDN0001', 'GDN0002', 'GDN0003', 'GDN0004', 'GDN0005', 'GDN0006', 'GDN0007', 'GDN0008', 'GDN0009', 'GDN0010', 'GDN0011', 'GDN0012', 'GDN0013', 'GDN0014', 'GDN0015', 'GDN0016', 'GDN0017', 'GDN0018', 'GDN0019', 'GDN0020', 'GDN0021', 'GDN0022', 'GDN0023', 'GDN0024', 'GDN0025', 'GDN0026', 'GDN0027', 'GDN0028', 'GDN0029', 'GDN0030'};
gender = {'Female', 'Female', 'Male', 'Female','Female','Female','Female','Male','Male','Male', 'Female', 'Male', 'Female','Female', 'Male', 'Female','Female','Male', 'Female','Female','Male', 'Female','Male', 'Female','Male', 'Female','Male','Male','Male','Male'};

% Fundamental Frequency at Resting Position
fundamental_frequency_gender = [0.120225 0.121424 0.095182 0.165967 0.140703 0.19972 0.107099 0.113163 0.127642 0.213669 0.3036 0.15211 0.214791 0.210433 0.112306 0.249472 0.178956 0.192393 0.240135 0.185054 0.295094 0.227227 0.129669 0.24572 0.0692 0.181502 0.182041 0.12256 0.161082 0.10576];

% dB Overtone at Resting Position
db_overtone = [-1.579378 -25.283933 -2.04023 -1.044529 -4.759265 -4.923702 -4.408474 -12.07235 -10.349847 -1.846589 -2.345343 -0.747955 -7.53673 -0.757966 -2.197296 -1.25706 -0.827158 -1.469303 -1.693414 -1.351698 -0.334193 -1.660997 -1.477652 -2.521721 -8.139457 -1.444331 -0.489974 -3.23775 -1.143744 -5.792813];
db_overtone = abs(db_overtone);
%% Visualization Funadametal Frequency for Gender Classification Analysis

% Create a scatter plot
figure;
scatter(1:length(subjects_gender), fundamental_frequency_gender, [], 'filled');

hold on;

% Separate data based on gender
female_indices = strcmp(gender, 'Female');
male_indices = strcmp(gender, 'Male');

% Plot female data points
scatter(find(female_indices), fundamental_frequency_gender(female_indices), 'r', 'filled');

% Plot male data points
scatter(find(male_indices), fundamental_frequency_gender(male_indices), 'y', 'filled');

% Add threshold lines
line([1, length(subjects_gender)], [0.165, 0.165], 'Color', 'b', 'LineStyle', '--'); 
line([1, length(subjects_gender)], [0.31, 0.31], 'Color', 'b', 'LineStyle', '--');

% Set x-axis labels
xticks(1:length(subjects_gender));
xticklabels(subjects_gender);
xlabel('Subject ID');

% Set y-axis label
ylabel('Fundamental Frequency');

% Add legend
legend('Threshold', 'Female', 'Male', 'Location', 'NorthEast');

% Title
title('Fundamental Frequency by Gender');

% Add grid
grid on;

%% Visualization dB Overtone for Gender Classification Analysis

% Create a scatter plot
figure;
scatter(1:length(subjects_gender), db_overtone, [], 'filled');

hold on;

% Separate data based on gender
female_indices = strcmp(gender, 'Female');
male_indices = strcmp(gender, 'Male');

% Plot female data points
scatter(find(female_indices), db_overtone(female_indices), 'r', 'filled');

% Plot male data points
scatter(find(male_indices), db_overtone(male_indices), 'y', 'filled');

% Set x-axis labels
xticks(1:length(subjects_gender));
xticklabels(subjects_gender);
xlabel('Subject ID');

% Set y-axis label
ylabel('dB Overtone');

% Add legend
legend('Female', 'Male', 'Location', 'NorthEast');

% Title
title('dB Overtone by Gender');

% Add grid
grid on;


