% Run Truss Tcl file in opensees
clear
close
clc

%% DEFINE INPTUTS
% Primary Inputs
analysis.model_id = 1;
analysis.gm_id = 4;
analysis.name = 'test';

% Secondary Inputs
analysis.type = 3;
analysis.max_disp = 1;
analysis.time_step = 0.01;
analysis.nonlinear = 2;

tic
%% Initial Setup
import tools.*
import display_model.model_plot

%% Load Analysis and Model parameters
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
ground_motion = gm_table(gm_table.id == analysis.gm_id,:);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
% Create Outputs Directory
output_dir = ['outputs/' model.name{1} '/' analysis.name];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Create Model Databases
[ node, element, story, joint, wall, hinge ] = fn_model_table( model, analysis );

%% Write TCL file
[ node ] = fn_build_model( output_dir, node, element, story, joint, wall, hinge, analysis );
fn_define_recorders( output_dir, analysis.type, node.id )
fn_define_loads( output_dir, analysis, model.damp_ratio, node, ground_motion )
% fn_eigen_analysis( output_dir, analysis.time_step, story.first_story_node )
fn_define_analysis( output_dir, analysis, node.id, ground_motion.eq_length, ground_motion.eq_dt )

%% Run Opensees
command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
system(command);

%% Plot Model in Matlab for verification
% model_plot('3D',1,0,1,1,1,0);

%% Save workspace data
save([output_dir filesep 'analysis_data'])

toc

