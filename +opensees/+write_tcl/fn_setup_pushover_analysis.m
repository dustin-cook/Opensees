function [ ] = fn_setup_pushover_analysis( output_dir, analysis, primary_nodes, story )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    
%% Write Pushover Analysis File
file_name = [output_dir filesep 'setup_pushover_analysis.tcl'];
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
tolerance = 5e-5;
fprintf(fileID,'test NormDispIncr %f 1000 \n',tolerance);

% Define Solution Algorithm
fprintf(fileID,'algorithm KrylovNewton \n');

% Define Each Load Step
control_node = primary_nodes(end);
if strcmp(analysis.pushover_direction,'x')
    control_dof = 1;
elseif strcmp(analysis.pushover_direction,'z')
    control_dof = 3;
end
max_displacement = analysis.pushover_drift*(story.y_start(end)+story.story_ht(end));
step_size = max_displacement / analysis.pushover_num_steps;
int_controller = ['DisplacementControl ' num2str(control_node) ' ' num2str(control_dof) ' ' num2str(step_size)]; 
fprintf(fileID,'integrator %s \n',int_controller);

% Define analysis type
fprintf(fileID,'analysis Static \n');

% Close File
fclose(fileID);


end

