% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup 
import tools.*

%% Define Model Parameters
% 2d linear model single frame, 1 bay, uniform members
num_bays = 3;
num_stories = 2;
story_ht_in = 120;
bay_width_in = 240;
foundation_fix = [1 1 1];
story_mass = 5;
E = 60000;
A = 288;
I = 13824;
story_force_k = 0;
story_weight_k = 0;
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

% Create Model Databases
[ node, element ] = fn_model_table( num_stories, num_bays, story_ht_in, bay_width_in, foundation_fix, story_mass, story_weight_k, story_force_k, A, E, I );

%% Write TCL file
fn_build_model( analysis.id, node, element )
fn_define_recorders( analysis.id, node.id, element.id )
fn_define_loads( analysis.id, damp_ratio, analysis.type, node )
fn_define_analysis( analysis, node.id )

%% Run Opensees
command = ['opensees ' 'Analysis' filesep analysis.id filesep 'run_analysis.tcl'];
system(command);

%% Load outputs and plot
% Nodal Displacement
node.disp_x = dlmread([output_local filesep 'nodal_disp_x.txt'],' ')';
node.disp_y = dlmread([output_local filesep 'nodal_disp_y.txt'],' ')';

% Nodal Reaction
node.reac_x = dlmread([output_local filesep 'nodal_reaction_x.txt'],' ')';
node.reac_y = dlmread([output_local filesep 'nodal_reaction_y.txt'],' ')';

% Element Forces
for i = 1:length(element.id)
ele_force = dlmread([output_local filesep 'element_' num2str(i) '_force.txt'],' ');
element.fx1{i} = ele_force(:,1)';
element.fy1{i} = ele_force(:,2)';
element.mz1{i} = ele_force(:,3)';
element.fx2{i} = ele_force(:,4)';
element.fy2{i} = ele_force(:,5)';
element.mz2{i} = ele_force(:,6)';
end

%% Display Results
if analysis.type == 1 % static load analysis
    disp(element)
elseif analysis.type == 2 % pushover analysis
    if num_bays ==0
        plot(abs(node.disp_x(end,:)),abs(node.reac_x(1,:)));
    else
        plot(abs(node.disp_x(end,:)),abs(sum(node.reac_x(1:(num_bays+1),:))));
    end
    grid on
elseif analysis.type == 3 % dynamic analysis
    plot(node.disp_x(end,:))
else
    error('Unkown Analysis Type')
end

toc
