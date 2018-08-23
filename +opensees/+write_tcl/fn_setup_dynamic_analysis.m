function [ ] = fn_setup_dynamic_analysis( output_dir, analysis )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Write Dynamic Analysis File
file_name = [output_dir filesep 'setup_dynamic_analysis.tcl'];
fileID = fopen(file_name,'w');

% Clear set up for this analysis
fprintf(fileID,'wipe \n');

% Build Model and Analysis Parameters
fprintf(fileID,'source %s/model.tcl \n', output_dir);
if analysis.run_eigen
    fprintf(fileID,'source %s/eigen.tcl \n', output_dir);
end
fprintf(fileID,'source %s/loads.tcl \n', output_dir);
fprintf(fileID,'source %s/recorders.tcl \n', output_dir);

% ANALYSIS DEFINITION
fprintf(fileID,'wipeAnalysis \n');

% Define Constraints
fprintf(fileID,'constraints Transformation \n');

% Define the DOF_numbered object
fprintf(fileID,'numberer RCM \n');

% Construct Linear Solver and linear SOE Objects
fprintf(fileID,'system BandGeneral \n');

% Test for Convergence
tolerance = 1e-6;
fprintf(fileID,'test NormDispIncr %f 1000 \n',tolerance);

% Define Solution Algorithm
fprintf(fileID,'algorithm NewtonLineSearch \n');

% Define Each Load Step
fprintf(fileID,'integrator Newmark 0.5 0.25 \n');

% Define analysis type
fprintf(fileID,'analysis Transient \n');

% Close File
fclose(fileID);

end
