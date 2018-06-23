% Combine DCR from load combos and plot
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '09DL';
additional_load_cases = {'11DL11LL'};

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
load([output_dir filesep 'ASCE_data']);

% For each Load case pull results
for i = 1:length(additional_load_cases)
    % Elements
    output_dir = ['outputs' filesep model.name{1} filesep additional_load_cases{i}];
    load_cases{i} = load([output_dir filesep 'ASCE_data']);
end

%% Calculate Envelopes
for i = 1:length(additional_load_cases)
    % DCR
    element.DCR_M = max([element.DCR_M,load_cases{i}.element.DCR_M],[],2);
    element.DCR_V = max([element.DCR_V,load_cases{i}.element.DCR_V],[],2);
    element.DCR_P = max([element.DCR_P,load_cases{i}.element.DCR_P],[],2);
    
    % EDP Profiles
    for j = 1:length(dirs_ran)
        story.(['max_accel_' dirs_ran(j)]) = max([story.(['max_accel_' dirs_ran(j)]),load_cases{i}.story.(['max_accel_' dirs_ran(j)])],[],2);
        story.(['max_disp_' dirs_ran(j)]) = max([story.(['max_disp_' dirs_ran(j)]),load_cases{i}.story.(['max_disp_' dirs_ran(j)])],[],2);
        story.(['max_drift_' dirs_ran(j)]) = max([story.(['max_drift_' dirs_ran(j)]),load_cases{i}.story.(['max_drift_' dirs_ran(j)])],[],2);
    end
end

% Total Envelope DCR
element.DCR_total = max([element.DCR_M,element.DCR_V,element.DCR_P],[],2);

%% Save Data
combo_dir = ['outputs' filesep model.name{1} filesep 'load_combo'];

% remove analysis data
clear analysis
clear output_dir

% Check to see if plot dir exists 
if ~exist(combo_dir,'dir')
  mkdir(combo_dir)  
end

% save data
save([combo_dir filesep 'ASCE_data'])
