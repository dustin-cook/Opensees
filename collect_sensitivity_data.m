clear all
close all
clc

% Collect results from model sensitivity study

% Define Model
analysis.model_id = 18;
analysis.proceedure = 'NDP';
analysis.id = 'baseline_beams_include';
analysis.model = 'ICBS_model_5ew_col_base';

%% Define outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study' ];
 
%% Define sensitivity study parameters
% model_name = {'ductility', 'ductility_no_var', 'strength', 'both', 'both_more'};
% model_name = {'ductility', 'strength', 'both'};
model_name = {'both'};
num_bays = 5; % just hardcoded for now

%% Collect baseline data
% Define model directories
baseline_dir = ['outputs' filesep analysis.model filesep analysis_name];
baseline_input_dir = [baseline_dir filesep 'asce_41_data'];

% Calculate input parameters
load([baseline_input_dir filesep 'element_analysis.mat'])
load([baseline_input_dir filesep 'joint_analysis.mat'])
load([baseline_input_dir filesep 'story_analysis.mat'])
[strength_cov, ductility_cov, scwb_ratio] = collect_model_inputs(element, story, joint);

% Load sensitivity summary results
summary_results = readtable([baseline_dir filesep 'IDA' filesep 'Fragility Data' filesep 'summary_outputs.csv']);

% Create Outputs Table
id = 1;
[outputs_table] = create_outputs_table(id, 'baseline', 1, strength_cov, ductility_cov, scwb_ratio, num_bays, summary_results);

%% Collect data for each sensitivity study
for m = 1:length(model_name)
    % Read all models
    model_dir = [outputs_dir filesep model_name{m} filesep 'model_files'];
    
    models = dir([model_dir filesep 'model_*']);
    for mdl = 1:length(models)
        % Define model directories
        this_model_dir = [model_dir filesep models(mdl).name];
        inputs_dir = [this_model_dir filesep 'asce_41_data'];
        load([inputs_dir filesep 'element_analysis.mat'])
        load([inputs_dir filesep 'joint_analysis.mat'])
        load([baseline_input_dir filesep 'story_analysis.mat'])

        % Calculate input parameters
        [strength_cov, ductility_cov, scwb_ratio] = collect_model_inputs(element, story, joint);

        % Load sensitivity summary results
        summary_results = readtable([this_model_dir filesep 'IDA' filesep 'Fragility Data' filesep 'summary_outputs.csv']);

        % Join tables
        id = id + 1;
        [outputs_table_row] = create_outputs_table(id, model_name{m}, mdl, strength_cov, ductility_cov, scwb_ratio, num_bays, summary_results);
        outputs_table = [outputs_table; outputs_table_row];
    end
end

% Write sensitivity study outputs file
writetable(outputs_table,[outputs_dir filesep 'summary_results.csv'])


% Functions 
function [strength_cov, ductility_cov, scwb_ratio] = collect_model_inputs(element, story, joint)
    first_story_columns = element(element.story == 1 & strcmp(element.type,'column'),:);
    strength_cov = std(first_story_columns.Mn_pos_1)/mean(first_story_columns.Mn_pos_1);
    ductility_cov = std(first_story_columns.b_hinge_1)/mean(first_story_columns.b_hinge_1);
    for s = 1:height(story)
        story.scwb(s) = mean(joint.col_bm_ratio(joint.story == s));
    end
    scwb_ratio = mean(story.scwb);
end

function [outputs_table_row] = create_outputs_table(id, group, variant, strength_cov, ductility_cov, scwb_ratio, num_bays, summary_results)
    inputs.id = id;
    summary_results.id = id;
    inputs.group = {group};
    inputs.variant = variant;
    inputs.strength_cov = strength_cov;
    inputs.ductility_cov = ductility_cov;
    inputs.scwb_ratio = scwb_ratio;
    inputs.num_bays = num_bays;
    inputs_table = struct2table(inputs);
    outputs_table_row = join(inputs_table, summary_results);
end

