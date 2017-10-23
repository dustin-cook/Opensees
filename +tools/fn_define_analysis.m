function [ ] = fn_define_analysis( analysis, nodes )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Define Parameters
if analysis.type == 1 % static load analysis
    int_controller = 'LoadControl 1';
    num_steps = 1;
    analysis_str_id = 'Static';
    time_step = 0;
elseif analysis.type == 2 % pushover analysis
    control_node = nodes(end);
    control_dof = 1;
    num_steps = 10;
    step_size = analysis.max_displ / num_steps;
    int_controller = ['DisplacementControl ' num2str(control_node) ' ' num2str(control_dof) ' ' num2str(step_size)]; 
    analysis_str_id = 'Static';
    time_step = 0;
elseif analysis.type == 3  % dynamic analysis
    gamma = 0.5; %trapezoidal
    beta = 0.25;
    int_controller = ['Newmark' ' ' num2str(gamma) ' ' num2str(beta)];
    analysis_str_id = 'Transient';
    num_steps = 3995;
    time_step = 0.1;
else
    error('Unkown Analysis Type')
end

%% Write Loads File
file_name = ['Analysis' filesep analysis.id filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

% Clear set up for this analysis
fprintf(fileID,'## Clear set up for this analysis \n');
fprintf(fileID,'wipe \n');

% Build Model and Analysis Parameters
fprintf(fileID,'## Build Model and Analysis Parameters \n');
fprintf(fileID,'source Analysis/%s/model.tcl \n', analysis.id);
fprintf(fileID,'source Analysis/%s/recorders.tcl \n', analysis.id);
fprintf(fileID,'source Analysis/%s/loads.tcl \n', analysis.id);

% ANALYSIS DEFINITION
fprintf(fileID,'## ANALYSIS DEFINITION \n');

% Define Constraints
fprintf(fileID,'# Define Constraints \n');
fprintf(fileID,'constraints Plain \n');

% Define the DOF_numbered object
fprintf(fileID,'# Define the DOF_numbered object \n');
fprintf(fileID,'numberer Plain \n');

% Construct Linear Solver and linear SOE Objects
fprintf(fileID,'# Construct Linear Solver and linear SOE Objects \n');
fprintf(fileID,'system BandGeneral \n');

% % Construct Convergence Test
% fprintf(fileID,'# Construct Convergence Test \n');
% fprintf(fileID,'test NormDispIncr 1.0e-6 6 \n');

% Define Solution Algorithm
fprintf(fileID,'# Define Solution Algorithm \n');
fprintf(fileID,'algorithm Linear \n');

% Define Each Load Step (displacement controlled)
fprintf(fileID,'# Define Each Load Step (displacement controlled) \n');
fprintf(fileID,'integrator %s \n',int_controller);

% Define analysis type
fprintf(fileID,'# Define analysis type \n');
fprintf(fileID,'analysis %s \n',analysis_str_id);

%% Run the Analysis
fprintf(fileID,'## Run the Analysis \n');
fprintf(fileID,'analyze %d %d \n',num_steps, time_step);
fprintf(fileID,'puts "Done!" \n');
fprintf(fileID,'wipe \n');

% Close File
fclose(fileID);


end

