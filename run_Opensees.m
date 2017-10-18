% Run Truss Tcl file in opensees
clear
close
clc

%% Initial Setup 
import tools.*

%% Define Model Parameters
% 2d linear model single frame, 1 bay, uniform members
story_ht_in = 120;
bay_width_in = 240;
foundation_fix = [1 1 0];
story_mass = 0;
E = 9999999999;
A = 9999999999;
I = 9999999999;
story_force_k = 10;
story_wieght_k = 0;

% Create Outputs Directory
analysis_id = 'test';
output_local = ['Analysis' filesep analysis_id];
if ~exist(output_local,'dir')
    mkdir(output_local);
end

%% Write TCL file
fn_build_model( story_ht_in, bay_width_in, foundation_fix, story_mass, E, A, I, analysis_id )
fn_define_recorders( analysis_id )
fn_define_loads( story_force_k, story_wieght_k, analysis_id )
fn_define_analysis( analysis_id )

%% Run Opensees
command = ['opensees ' 'Analysis' filesep analysis_id filesep 'run_analysis.tcl'];
system(command);

%% Load outputs and plot
nodal_disp = readtable([output_local filesep 'nodal_disp_x.txt'],'Delimiter',' ','ReadVariableNames',false);
element_1_force = readtable([output_local filesep 'element_1_force.txt'],'Delimiter',' ','ReadVariableNames',false);
element_2_force = readtable([output_local filesep 'element_2_force.txt'],'Delimiter',' ','ReadVariableNames',false);
element_3_force = readtable([output_local filesep 'element_3_force.txt'],'Delimiter',' ','ReadVariableNames',false);
