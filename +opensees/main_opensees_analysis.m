function [ ] = main_opensees_analysis( analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
import opensees.*

% Load inital Model info
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
% Create Outputs Directory
output_dir = ['outputs/' model.name{1} '/' analysis.name];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end 

%% Load Model
node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
story = readtable([output_dir filesep 'story.csv'],'ReadVariableNames',true);
joint = readtable([output_dir filesep 'joint.csv'],'ReadVariableNames',true);
truss = readtable([output_dir filesep 'truss.csv'],'ReadVariableNames',true);
if analysis.nonlinear == 0
    hinge = [];
else
    hinge = readtable([output_dir filesep 'hinge.csv'],'ReadVariableNames',true);
end

%% Write TCL file
[ node, ground_motion ] = main_write_tcl( model, output_dir, node, element, story, joint, hinge, analysis, truss );

%% Run Opensees
if analysis.run_opensees
    tic
    main_run_opensees( output_dir )
    toc
end

%% Postprocess OS data
main_post_process_opensees( analysis, model, story, node, element, ground_motion, output_dir )

end

