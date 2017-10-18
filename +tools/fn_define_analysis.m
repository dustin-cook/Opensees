function [ ] = fn_define_analysis( analysis_id )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Write Loads File
file_name = ['Analysis' filesep analysis_id filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

%% Clear set up for this analysis
fprintf(fileID,'## Clear set up for this analysis \n');
fprintf(fileID,'wipe \n');

%% Build Model and Analysis Parameters
fprintf(fileID,'## Build Model and Analysis Parameters \n');
fprintf(fileID,'source Analysis/%s/model.tcl \n', analysis_id);
fprintf(fileID,'source Analysis/%s/recorders.tcl \n', analysis_id);
fprintf(fileID,'source Analysis/%s/loads.tcl \n', analysis_id);

%% ANALYSIS DEFINITION
fprintf(fileID,'## ANALYSIS DEFINITION \n');

% Define Constraints
fprintf(fileID,'# Define Constraints \n');
fprintf(fileID,'constraints Transformation \n');

% Define the DOF_numbered object
fprintf(fileID,'# Define the DOF_numbered object \n');
fprintf(fileID,'numberer RCM \n');

% Construct Linear Solver and linear SOE Objects
fprintf(fileID,'# Construct Linear Solver and linear SOE Objects \n');
fprintf(fileID,'system BandGeneral \n');

% Construct Convergence Test
fprintf(fileID,'# Construct Convergence Test \n');
fprintf(fileID,'test NormDispIncr 1.0e-6 6 \n');

% Define Solution ALgorithm
fprintf(fileID,'# Define Solution ALgorithm \n');
fprintf(fileID,'algorithm Newton \n');

% Define Each Load Step (displacement controlled)
fprintf(fileID,'# Define Each Load Step (displacement controlled) \n');
fprintf(fileID,'integrator LoadControl 1 \n');

% Define analysis type
fprintf(fileID,'# Define analysis type \n');
fprintf(fileID,'analysis Static \n');

%% Run the Analysis
fprintf(fileID,'## Run the Analysis \n');
fprintf(fileID,'analyze 1 \n');

% Close File
fclose(fileID);

end

