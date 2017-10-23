% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup 
import tools.*

%% Define Model Parameters
% 2d linear model single frame, 1 bay, uniform members
num_bay = 1;
story_ht_in = 120;
bay_width_in = 240;
foundation_fix = [1 1 1];
story_mass = 0;
E = 60000;
A = 288;
I = 13824;
story_force_k = 0;
story_wieght_k = 0;
damp_ratio = 0.05;

%% Define Analysis parameters
analysis.id = 'test';
analysis.type = 3; % 1 = static force analysis, 2 = pushover analysis
analysis.max_displ = 1; % only for pushover analsys

% Create Outputs Directory
output_local = ['Analysis' filesep analysis.id];
if ~exist(output_local,'dir')
    mkdir(output_local);
end

%% Write TCL file
fn_build_model( story_ht_in, bay_width_in, foundation_fix, story_mass, E, A, I, analysis.id )
fn_define_recorders( analysis.id )
fn_define_loads( story_force_k, story_wieght_k, analysis.id, damp_ratio, analysis.type )
fn_define_analysis( analysis )

%% Run Opensees
command = ['opensees ' 'Analysis' filesep analysis.id filesep 'run_analysis.tcl'];
system(command);

%% Load outputs and plot
% Nodal Displacement
nodal_disp_x = dlmread([output_local filesep 'nodal_disp_x.txt'],' ')';
nodal_disp_y = dlmread([output_local filesep 'nodal_disp_y.txt'],' ')';
node = (1:length(nodal_disp_x(:,1)))';
nodal_disp = table(node, nodal_disp_x, nodal_disp_y, 'VariableNames', {'node','disp_x','disp_y'});

% Nodal Reaction
nodal_reac_x = dlmread([output_local filesep 'nodal_reaction_x.txt'],' ')';
nodal_reac_y = dlmread([output_local filesep 'nodal_reaction_y.txt'],' ')';
node = (1:length(nodal_reac_x(:,1)))';
nodal_reac = table(node, nodal_reac_x, nodal_reac_y, 'VariableNames', {'node','reac_x','reac_y'});

% Element Forces
ele_force_1 = dlmread([output_local filesep 'element_1_force.txt'],' ');
ele_force_2 = dlmread([output_local filesep 'element_2_force.txt'],' ');
ele_force_3 = dlmread([output_local filesep 'element_3_force.txt'],' ');

fx1 = [ele_force_1(:,1),ele_force_2(:,1),ele_force_3(:,1)]';
fy1 = [ele_force_1(:,2),ele_force_2(:,2),ele_force_3(:,2)]';
mz1 = [ele_force_1(:,3),ele_force_2(:,3),ele_force_3(:,3)]';
fx2 = [ele_force_1(:,4),ele_force_2(:,4),ele_force_3(:,4)]';
fy2 = [ele_force_1(:,5),ele_force_2(:,5),ele_force_3(:,5)]';
mz2 = [ele_force_1(:,6),ele_force_2(:,6),ele_force_3(:,6)]';

element = (1:3)';
element_force = table(element,fx1,fy1,mz1,fx2,fy2,mz2, 'VariableNames', {'element','fx1','fy1','mz1','fx2','fy2','mz2'});

% Display Results
%% Define Parameters
if analysis.type == 1 % static load analysis
    disp(element_force)
elseif analysis.type == 2 % pushover analysis
    plot(abs(nodal_disp.disp_x(3,:)),abs(sum(nodal_reac.reac_x)));
    grid on
else
    error('Unkown Analysis Type')
end


toc
