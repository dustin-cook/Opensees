% Run Truss Tcl file in opensees
clear
close
clc

%% DEFINE INPTUTS
% Primary Inputs
analysis.model_id = 6;
analysis.gm_seq_id = 10;
analysis.name = 'scale_4';

% Secondary Inputs
analysis.type = 3;
analysis.max_disp = 7;
analysis.time_step = 0.005;
analysis.nonlinear = 2;
analysis.dead_load = 1.0;
analysis.live_load = 0;
analysis.accidental_torsion = 0;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.01;
analysis.hinge_stiff_mod = 10;

tic
%% Initial Setup
import tools.*
import display_model.model_plot

%% Load Analysis and Model parameters
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
% Create Outputs Directory
output_dir = ['outputs/' model.name{1} '/' analysis.name];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Create Model Databases
[ node, element, story, joint, hinge, truss ] = fn_model_table( model, analysis );
% Save element and node databases
writetable(node,[output_dir filesep 'node.csv'])
writetable(element,[output_dir filesep 'element.csv'])

%% Write TCL file
if strcmp(model.dimension,'2D')
    [ node ] = fn_build_model_2D( output_dir, node, element, story, joint, hinge, analysis, truss );
elseif strcmp(model.dimension,'3D')
    [ node ] = fn_build_model_3D( output_dir, node, element, story, joint, hinge, analysis );
else
    error('Number of Dimensions Not Recognized')
end
fn_define_recorders( output_dir, model.dimension, node.id', element, hinge )
[ground_motion] = fn_define_loads( output_dir, analysis, node, model.dimension, length(story.id), element.id');
% fn_eigen_analysis( output_dir, analysis.time_step, story.first_story_node, length(story.id))
fn_define_analysis( output_dir, analysis, node.id, ground_motion )

%% Run Opensees
command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
system(command);

%% Plot Model in Matlab for verification
% model_plot('3D',1,0,1,1,1,0);

%% Save workspace data
save([output_dir filesep 'analysis_data'])

toc

