function [ model ] = main_eigen_analysis( model, analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import opensees.*
import opensees.write_tcl.*

% Create Write Directory
% TCL or Opensees does not like filesep command on windows, therefore must manually define forward slash seperators
write_dir_opensees = [strrep(analysis.out_dir,'\','/') '/eigen_analysis']; 
fn_make_directory( write_dir_opensees )

% Define Read Directories
read_dir_model = [analysis.out_dir filesep 'model_data'];

% Load Model Data
node = readtable([read_dir_model filesep 'node.csv'],'ReadVariableNames',true);
element = readtable([read_dir_model filesep 'element.csv'],'ReadVariableNames',true);
story = readtable([read_dir_model filesep 'story.csv'],'ReadVariableNames',true);
joint = readtable([read_dir_model filesep 'joint.csv'],'ReadVariableNames',true);
hinge = readtable([read_dir_model filesep 'hinge.csv'],'ReadVariableNames',true);

%% Run Opensees for Eigen Analysis
% Write TCL files
[ ~ ] = fn_define_model( write_dir_opensees, node, element, joint, hinge, analysis, model.dimension, story, [], model );
primary_nodes = node.id(node.primary_story == 1 & node.story > 0);
fn_eigen_analysis( write_dir_opensees, primary_nodes', max(story.id), analysis, model.dimension)
fn_define_eignen_run_script( write_dir_opensees )

% Run Opensees
fprintf('Running Opensees Eigen Analysis\n')
main_run_opensees( write_dir_opensees, analysis )

% Postprocess OS data
periods = dlmread([write_dir_opensees filesep 'period.txt']);
model.T1_x = periods(1);
writetable(model,[read_dir_model filesep 'model.csv'])

end

