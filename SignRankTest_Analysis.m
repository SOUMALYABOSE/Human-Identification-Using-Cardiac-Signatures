%% Important Copyright Information with regards to this scirpt:
% Script has been added under the CME-Master Thesis titled "Human Identification With Gender Knowledge Analysis From Cardiac Signatures and Health-related Information Extraction From Radar Based Biomarkers" 
% Conducted @ Home Automation Lab, Institute FAPS, FAU Erlangen-Nuremberg by Mr. Soumalya Bose
% Advisors: Mr. Jochen Bauer (Institute FAPS, FAU), Dr. -med. Tobias Steigleder (PallMeT, Univ. Hospital Erlangen, FAU), Professor Dr.-Ing. Jörg Franke (Institute FAPS, FAU), Professor Dr.-Ing. Georg Fischer (Institute for Electronics Engineering, FAU)
% Copyright (C) 2023-24  Mr. Soumalya Bose
%% Clearing of Workspace
clear;
clc;

%% Declaring all arrays for Wilcoxon Rank-Sum test (WRST)
age = [24 24 24 38 38 38 25 25 25 28 28 28 27 27 27 49 49 49 24 24 24 24 24 24 40 40 40 24 24 24 21 21 21 24 24 24 49 49 49 35 35 35 27 27 27 42 42 42 49 49 49 28 28 28 27 27 27 23 23 23 24 24 24 61 61 61 27 27 27 21 21 21 26 26 26 31 31 31 24 24 24 29 29 29 25 25 25 25 25 25];

height = [166 166 166 161 161 161 187 187 187 178 178 178 173 173 173 172 172 172 187 187 187 182 182 182 184 184 184 186 186 186 165 165 165 193 193 193 173 173 173 153 153 153 182 182 182 165 165 165 167 167 167 165 165 165 166 166 166 172 172 172 187 187 187 178 178 178 186 186 186 165 165 165 183 183 183 160 160 160 187 187 187 190 190 190 186 186 186 172 172 172];
height = height ./100; % In Metres

weight = [63 63 63 49 49 49 82 82 82 59 59 59 93 93 93 63 63 63 80 80 80 77 77 77 74 74 74 78 78 78 55 55 55 86 86 86 62 62 62 44 44 44 78 78 78 61 61 61 85 85 85 57 57 57 59 59 59 68 68 68 85 85 85 90 90 90 69 69 69 65 65 65 82 82 82 51 51 51 83 83 83 94 94 94 82 82 82 93 93 93]; 

bmi = [28.9 28.9 28.9 18.7 18.7 18.7 23.4 23.4 23.4 18.6 18.6 18.6 31.1 31.1 31.1 21.3 21.3 21.3 22.9 22.9 22.9 23.2 23.2 23.2 21.9 21.9 21.9 22.5 22.5 22.5 20.2 20.2 20.2 23.1 23.1 23.1 20.7 20.7 20.7 18.8 18.8 18.8 23.5 23.5 23.5 22.4 22.4 22.4 30.5 30.5 30.5 20.9 20.9 20.9 21.4 21.4 21.4 23 23 23 24.3 24.3 24.3 28.4 28.4 28.4 19.9 19.9 19.9 23.9 23.9 23.9 24.5 24.5 24.5 19.7 19.7 19.7 23.7 23.7 23.7 26 26 26 23.7 23.7 23.7 31.4 31.4 31.4];

fundamental_frequency = [0.120225 0.050534 0 0.121424 0.250635 0 0.095182 0.098696 0 0.165967 0.106591 0.214479 0.140703 0.124331 0.185615 0.19972 0.088932 0.145233 0.107099 0.106166 0.102 0.113163 0.106675 0.121435 0.127642 0.131901 0.10552 0.213669 0.126719 0.14724 0.3036 0.151629 0.28308 0.15211 0.142773 0.137516 0.214791 0.12482 0.192872 0.210433 0.178754 0.135049 0.112306 0 0 0.249472 0.215832 0.167781 0.178956 0.098446 0.103734 0.192393 0.20162 0.115263 0.240135 0.260406 0.12098 0.185054 0.160868 0.186417 0.295094 0.253133 0.237866 0.227227 0.200046 0.128912 0.129669 0.135504 0.14744 0.24572 0 0 0.0692 0.090373 0.156949 0.181502 0 0 0.182041 0.078338 0.18488 0.12256 0.102753 0.106567 0.161082 0.151749 0.168683 0.10576 0.129469 0.147169];
fundamental_frequency = 2000 ./ fundamental_frequency; % In Samples/seconds
fundamental_frequency = fundamental_frequency ./ 1000;  % In Kilo Samples/seconds
fundamental_frequency(fundamental_frequency == Inf) = 0; % Replacing Infinite Value Elements With Zero Value

second_breathing_harmonic = [0.213957 0.225881 0 0.223347 0.505297 0 0.217843 0.20149 0 0.442731 0.237072 0.496937 0.23111 0.231769 0.377859 0.410852 0.20236 0.201844 0.245698 0.215829 0.289933 0.260276 0.202521 0.205848 0.206344 0.202729 0.211082 0.463133 0.217819 0.226565 0.602535 0.44161 0.649176 0.444136 0.219315 0.204082 0.42584 0.229828 0.413223 0.457332 0.420946 0.210553 0.206647 0 0 0.401278 0.414576 0.467475 0.493991 0.205592 0.268237 0.422973 0.412407 0.209589 0.452661 0.615412 0.201207 0.406155 0.444805 0.417384 0.622226 0.601955 0.410992 0.451686 0.452106 0.208164 0.218016 0.225715 0.215994 0.42426 0 0 0.277181 0.200345 0.41563 0.408027 0 0 0.444356 0.214012 0.453325 0.210819 0.21155 0.202354 0.308517 0.469644 0.471299 0.279285 0.242 0.230735];
second_breathing_harmonic = 2000 ./ second_breathing_harmonic; % In Samples/seconds
second_breathing_harmonic = second_breathing_harmonic ./ 1000;  % In Kilo Samples/seconds
second_breathing_harmonic(second_breathing_harmonic == Inf) = 0; % Replacing Infinite Value Elements With Zero Value

first_interHarmonic_distance = [0.093732 0.175347 0 0.101923 0.254662 0 0.122661 0.102793 0 0.276764 0.135996 0.282458 0.090406 0.107438 0.19224 0.211132 0.113427 0.056611 0.138599 0.109663 0.187933 0.147113 0.095847 0.084413 0.078702 0.070827 0.105562 0.249465 0.0911 0.079325 0.298935 0.28998 0.366096 0.292026 0.076542 0.066566 0.211048 0.105008 0.220351 0.246899 0.242192 0.075504 0.094341 0 0 0.151805 0.198743 0.299694 0.315035 0.107146 0.164503 0.23058 0.210787 0.094326 0.212526 0.355006 0.080227 0.221101 0.283936 0.230967 0.327132 0.348822 0.173125 0.224459 0.25206 0.079252 0.088347 0.090211 0.068554 0.17854 0 0 0.207981 0.109973 0.258681 0.226524 0 0 0.262315 0.135674 0.268445 0.088259 0.108797 0.095787 0.147435 0.317895 0.302616 0.173525 0.112631 0.083566];
first_interHarmonic_distance = 2000 ./ first_interHarmonic_distance; % In Samples/seconds
first_interHarmonic_distance = first_interHarmonic_distance ./ 1000;  % In Kilo Samples/seconds
first_interHarmonic_distance(first_interHarmonic_distance == Inf) = 0; % Replacing Infinite Value Elements With Zero Value

%% Perform WRST b/w Fundamental Frequency & Age
[p, h1, stats] = ranksum(fundamental_frequency, age);

% Display p-value and test statistic
disp(['p-value for Fund. Freq. and Age: ', num2str(p)]);
disp(['Test statistic for Fund. Freq. and Age: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([fundamental_frequency', age'], 'Labels', {'Fundamental Frequency (kSa/s)', 'Age'});
title('Boxplot of Fundamental Frequency and Age');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(fundamental_frequency), mean(age)]);
hold on;
errorbar([1, 2], [mean(fundamental_frequency), mean(age)], [std(fundamental_frequency), std(age)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Fundamental Frequency (kSa/s)', 'Age'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Fundamental Frequency & Height
[p, h2, stats] = ranksum(fundamental_frequency, height);

% Display p-value and test statistic
disp(['p-value for Fund. Freq. and Height: ', num2str(p)]);
disp(['Test statistic for Fund. Freq. and Height: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([fundamental_frequency', height'], 'Labels', {'Fundamental Frequency (kSa/s)', 'Height (m)'});
title('Boxplot of Fundamental Frequency and Height');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(fundamental_frequency), mean(height)]);
hold on;
errorbar([1, 2], [mean(fundamental_frequency), mean(height)], [std(fundamental_frequency), std(height)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Fundamental Frequency (kSa/s)', 'Height (m)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Fundamental Frequency & Weight
[p, h3, stats] = ranksum(fundamental_frequency, weight);

% Display p-value and test statistic
disp(['p-value for Fund. Freq. and Weight: ', num2str(p)]);
disp(['Test statistic for Fund. Freq. and Weight: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([fundamental_frequency', weight'], 'Labels', {'Fundamental Frequency (kSa/s)', 'Weight (kg)'});
title('Boxplot of Fundamental Frequency and Weight');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(fundamental_frequency), mean(weight)]);
hold on;
errorbar([1, 2], [mean(fundamental_frequency), mean(weight)], [std(fundamental_frequency), std(weight)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Fundamental Frequency (kSa/s)', 'Weight (kg)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Fundamental Frequency & BMI
[p, h4, stats] = ranksum(fundamental_frequency, bmi);

% Display p-value and test statistic
disp(['p-value for Fund. Freq. and BMI: ', num2str(p)]);
disp(['Test statistic for Fund. Freq. and BMI: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([fundamental_frequency', bmi'], 'Labels', {'Fundamental Frequency (kSa/s)', 'BMI (kg/m^2)'});
title('Boxplot of Fundamental Frequency and BMI');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(fundamental_frequency), mean(bmi)]);
hold on;
errorbar([1, 2], [mean(fundamental_frequency), mean(bmi)], [std(fundamental_frequency), std(bmi)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Fundamental Frequency (kSa/s)', 'BMI (kg/m^2)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Second Breathing Harmonic & Age
[p, h5, stats] = ranksum(second_breathing_harmonic, age);

% Display p-value and test statistic
disp(['p-value for Second Breathing Harmonic and Age: ', num2str(p)]);
disp(['Test statistic for Second Breathing Harmonic and Age: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([second_breathing_harmonic', age'], 'Labels', {'Second Breathing Harmonic (kSa/s)', 'Age'});
title('Boxplot of Second Breathing Harmonic and Age');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(second_breathing_harmonic), mean(age)]);
hold on;
errorbar([1, 2], [mean(second_breathing_harmonic), mean(age)], [std(second_breathing_harmonic), std(age)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Second Breathing Harmonic (kSa/s)', 'Age'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Second Breathing Harmonic & Height
[p, h6, stats] = ranksum(second_breathing_harmonic, height);

% Display p-value and test statistic
disp(['p-value for Second Breathing Harmonic and Height: ', num2str(p)]);
disp(['Test statistic for Second Breathing Harmonic and Height: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([second_breathing_harmonic', height'], 'Labels', {'Second Breathing Harmonic (kSa/s)', 'Height (m)'});
title('Boxplot of Second Breathing Harmonic and Height');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(second_breathing_harmonic), mean(height)]);
hold on;
errorbar([1, 2], [mean(second_breathing_harmonic), mean(height)], [std(second_breathing_harmonic), std(height)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Second Breathing Harmonic (kSa/s)', 'Height (m)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Second Breathing Harmonic & Weight
[p, h7, stats] = ranksum(second_breathing_harmonic, weight);

% Display p-value and test statistic
disp(['p-value for Second Breathing Harmonic and Weight: ', num2str(p)]);
disp(['Test statistic for Second Breathing Harmonic and Weight: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([second_breathing_harmonic', weight'], 'Labels', {'Second Breathing Harmonic (kSa/s)', 'Weight (kg)'});
title('Boxplot of Second Breathing Harmonic and Weight');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(second_breathing_harmonic), mean(weight)]);
hold on;
errorbar([1, 2], [mean(second_breathing_harmonic), mean(weight)], [std(second_breathing_harmonic), std(weight)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Second Breathing Harmonic (kSa/s)', 'Weight (kg)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w Second Breathing Harmonic & BMI
[p, h8, stats] = ranksum(second_breathing_harmonic, bmi);

% Display p-value and test statistic
disp(['p-value for Second Breathing Harmonic and BMI: ', num2str(p)]);
disp(['Test statistic for Second Breathing Harmonic and BMI: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([second_breathing_harmonic', bmi'], 'Labels', {'Second Breathing Harmonic (kSa/s)', 'BMI (kg/m^2)'});
title('Boxplot of Second Breathing Harmonic and BMI');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(second_breathing_harmonic), mean(bmi)]);
hold on;
errorbar([1, 2], [mean(second_breathing_harmonic), mean(bmi)], [std(second_breathing_harmonic), std(bmi)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'Second Breathing Harmonic (kSa/s)', 'BMI (kg/m^2)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w First Inter-Harmonic Distance & Age
[p, h9, stats] = ranksum(first_interHarmonic_distance, age);

% Display p-value and test statistic
disp(['p-value for First Inter-Harmonic Distance and Age: ', num2str(p)]);
disp(['Test statistic for First Inter-Harmonic Distance and Age: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([first_interHarmonic_distance', age'], 'Labels', {'First Inter-Harmonic Distance (kSa/s)', 'Age'});
title('Boxplot of First Inter-Harmonic Distance and Age');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(first_interHarmonic_distance), mean(age)]);
hold on;
errorbar([1, 2], [mean(first_interHarmonic_distance), mean(age)], [std(first_interHarmonic_distance), std(age)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'First Inter-Harmonic Distance (kSa/s)', 'Age'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w First Inter-Harmonic Distance & Height
[p, h10, stats] = ranksum(first_interHarmonic_distance, height);

% Display p-value and test statistic
disp(['p-value for First Inter-Harmonic Distance and Height: ', num2str(p)]);
disp(['Test statistic for First Inter-Harmonic Distance and Height: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([first_interHarmonic_distance', height'], 'Labels', {'First Inter-Harmonic Distance (kSa/s)', 'Height (m)'});
title('Boxplot of First Inter-Harmonic Distance and Height');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(first_interHarmonic_distance), mean(height)]);
hold on;
errorbar([1, 2], [mean(first_interHarmonic_distance), mean(height)], [std(first_interHarmonic_distance), std(height)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'First Inter-Harmonic Distance (kSa/s)', 'Height (m)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w First Inter-Harmonic Distance & Weight
[p, h11, stats] = ranksum(first_interHarmonic_distance, weight);

% Display p-value and test statistic
disp(['p-value for First Inter-Harmonic Distance and Weight: ', num2str(p)]);
disp(['Test statistic for First Inter-Harmonic Distance and Weight: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([first_interHarmonic_distance', weight'], 'Labels', {'First Inter-Harmonic Distance (kSa/s)', 'Weight (kg)'});
title('Boxplot of First Inter-Harmonic Distance and Weight');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(first_interHarmonic_distance), mean(weight)]);
hold on;
errorbar([1, 2], [mean(first_interHarmonic_distance), mean(weight)], [std(first_interHarmonic_distance), std(weight)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'First Inter-Harmonic Distance (kSa/s)', 'Weight (kg)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');

%% Perform WRST b/w First Inter-Harmonic Distance & BMI
[p, h12, stats] = ranksum(first_interHarmonic_distance, bmi);

% Display p-value and test statistic
disp(['p-value for First Inter-Harmonic Distance and BMI: ', num2str(p)]);
disp(['Test statistic for First Inter-Harmonic Distance and BMI: ', num2str(stats.ranksum)]);

% Plotting
figure;
subplot(2, 1, 1);
boxplot([first_interHarmonic_distance', bmi'], 'Labels', {'First Inter-Harmonic Distance (kSa/s)', 'BMI (kg/m^2)'});
title('Boxplot of First Inter-Harmonic Distance and BMI');
ylabel('Median Comparison');

subplot(2, 1, 2);
bar([1, 2], [mean(first_interHarmonic_distance), mean(bmi)]);
hold on;
errorbar([1, 2], [mean(first_interHarmonic_distance), mean(bmi)], [std(first_interHarmonic_distance), std(bmi)], 'k.', 'LineWidth', 1);
set(gca, 'XTick', [1, 2], 'XTickLabels', {'First Inter-Harmonic Distance (kSa/s)', 'BMI (kg/m^2)'});
title('Mean Comparison');
ylabel('Mean');
legend('Mean', 'Standard Deviation', 'Location', 'northwest');


