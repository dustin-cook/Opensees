function [ ] = main_opensees_analysis( analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
import opensees.*

%% Load data
% Load inital Model info
if analysis.model_type == 1 % SDOF
    model_table = readtable(['inputs' filesep 'sdof_models.csv'],'ReadVariableNames',true);
    model = model_table(model_table.id == analysis.model_id,:); 
elseif analysis.model_type == 2 % MDOF
    model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
    model = model_table(model_table.id == analysis.model_id,:);
end

% Load Model
node = readtable(['outputs/' model.name{1} filesep 'model data' filesep 'node.csv'],'ReadVariableNames',true);
element = readtable(['outputs/' model.name{1} filesep 'model data' filesep 'element.csv'],'ReadVariableNames',true);
story = readtable(['outputs/' model.name{1} filesep 'model data' filesep 'story.csv'],'ReadVariableNames',true);
joint = readtable(['outputs/' model.name{1} filesep 'model data' filesep 'joint.csv'],'ReadVariableNames',true);
hinge = readtable(['outputs/' model.name{1} filesep 'model data' filesep 'hinge.csv'],'ReadVariableNames',true);

%% Create Outputs Directory
output_dir = ['outputs/' model.name{1} '/' analysis.name];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end 

%% Write TCL files
[ node, ground_motion ] = main_write_tcl( model, output_dir, node, element, story, joint, hinge, analysis );

%% Run Opensees
if analysis.run_opensees
    tic
    main_run_opensees( output_dir )
    toc
end

%% Postprocess OS data
main_post_process_opensees( analysis, model, story, node, element, hinge, ground_motion, output_dir )

end

