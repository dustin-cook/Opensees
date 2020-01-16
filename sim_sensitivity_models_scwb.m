clear all
close all
clc

% Creates variations on pre-defined models 
% Only for 2D models
% Only modifies the first story columns

import asce_41.*
import plotting_tools.*
import opensees.write_tcl.*

rng(111389)

% Define Model
analysis.model_id = 19;
analysis.proceedure = 'NDP';
analysis.id = 'baseline_fix';
analysis.model = 'ICBS_model_5ew_col_base_no_mods';

% Secondary options
analysis.type = 1;
analysis.nonlinear = 1;
analysis.stories_nonlinear = 6; % Default to all modeling all stories as nonlinear when doing NDP
analysis.elastic_beams = 1;
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;
analysis.suppress_outputs = 1;
analysis.joint_explicit = 0;
analysis.joint_model = 1; % 1 = beam/column elements, 2 = joint 3D

%% Define inputs and outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
inputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'asce_41_data'];
os_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'opensees_data'];
os_model_csv_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'model_data'];

%% Define Element Properties
% Load Data
ele_prop_table = readtable(['inputs' filesep 'element.csv']);
load([inputs_dir filesep 'element_analysis.mat'])
load([inputs_dir filesep 'joint_analysis.mat'])
load([inputs_dir filesep 'hinge_analysis.mat'])
load([inputs_dir filesep 'node_analysis.mat'])
load([inputs_dir filesep 'story_analysis.mat'])
    
%% Define sensitivity study parameters
model_name = 'scwb';

%% Create Models
% Create outputs dir
outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study' filesep model_name];
new_models_dir = [outputs_dir filesep 'model_files'];
mkdir(outputs_dir)
mkdir(new_models_dir)

%% Get SCWB profile
for s = 1:height(story)
    story.scwb(s) = mean(joint.col_bm_ratio(joint.story == s));
end

% Modify beam strengths to create new SCWB ratio's
bm_str = [2, 1.5, 1, 0.75, 0.5];
element2save = element;
for b = 1:length(bm_str)
    bm_filt = strcmp(element.type,'beam');
    for s = 1:2
        element.(['Mn_pos_' num2str(s)])(bm_filt) = bm_str(b)*element2save.(['Mn_pos_' num2str(s)])(bm_filt);
        element.(['Mn_neg_' num2str(s)])(bm_filt) = bm_str(b)*element2save.(['Mn_neg_' num2str(s)])(bm_filt);
        element.(['Mp_pos_' num2str(s)])(bm_filt) = bm_str(b)*element2save.(['Mp_pos_' num2str(s)])(bm_filt);
        element.(['Mp_neg_' num2str(s)])(bm_filt) = bm_str(b)*element2save.(['Mp_neg_' num2str(s)])(bm_filt);
    end
    
    %% Create plots and tcl script for each model
    model_dir = [new_models_dir filesep 'model_' num2str(b)];
    mkdir(model_dir)

    % Copy Model Files
    mkdir([model_dir filesep 'asce_41_data']);
    mkdir([model_dir filesep 'opensees_data']);
    mkdir([model_dir filesep 'model_data']);
    copyfile([os_model_csv_dir filesep 'node.csv'], [model_dir filesep 'model_data'])
    copyfile([os_model_csv_dir filesep 'story.csv'], [model_dir filesep 'model_data'])
    copyfile([os_model_csv_dir filesep 'hinge.csv'], [model_dir filesep 'model_data'])
    
    % Re-calculate SCWB ratio
    for j = 1:height(joint)
        [ jnt ] = fn_joint_capacity( joint(j,:), element, ele_prop_table, [] );
        joint(j,:) = jnt;
    end
    for s = 1:height(story)
        story.scwb(s) = mean(joint.col_bm_ratio(joint.story == s));
    end
    mean(story.scwb)
    
    % Write Element Analysis Data file with new a and b values
    save([model_dir filesep 'asce_41_data' filesep 'story_analysis.mat'],'story')
    save([model_dir filesep 'asce_41_data' filesep 'element_analysis.mat'],'element')
    save([model_dir filesep 'asce_41_data' filesep 'joint_analysis.mat'],'joint')

    % Write model.tcl file with new values
    [ ~ ] = fn_define_model( [model_dir filesep 'opensees_data'], node, element, joint, hinge, analysis, '2D', story, [model_dir filesep 'asce_41_data'] );
end



