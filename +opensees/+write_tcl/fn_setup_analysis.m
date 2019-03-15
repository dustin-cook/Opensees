function [ ] = fn_setup_analysis( write_dir, model_dir, analysis, first_story_node, story )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import opensees.write_tcl.*

%% Write Analysis Setup File
file_name = [write_dir filesep 'setup_analysis.tcl'];
fileID = fopen(file_name,'w');

% Clear set up for this analysis
fprintf(fileID,'wipe \n');

% Build Model and Analysis Parameters
fprintf(fileID,'source %s/model.tcl \n', model_dir);
if analysis.run_eigen
    fprintf(fileID,'source %s/eigen.tcl \n', write_dir);
end
fprintf(fileID,'source %s/loads.tcl \n', write_dir);
fprintf(fileID,'source %s/recorders.tcl \n', write_dir);

% ANALYSIS DEFINITION
fprintf(fileID,'wipeAnalysis \n');

% Define Constraints
fprintf(fileID,'constraints Transformation \n');

% Define the DOF_numbered object
fprintf(fileID,'numberer RCM \n');

% Construct Linear Solver and linear SOE Objects
if analysis.opensees_SP
    fprintf(fileID,'system Mumps \n'); % Use Mumps for OpenseesSP
%     fprintf(fileID,'system BandGeneral \n');
else
    fprintf(fileID,'system BandGeneral \n');
end

% Convergence test
tolerance = 1e-6;
fprintf(fileID,'test NormDispIncr %f 100 \n',tolerance);

% Define analysis type
fprintf(fileID,'algorithm KrylovNewton \n');

if analysis.type == 1 % Dynamic Analysis
    % Define Each Load Step
%     fprintf(fileID,'integrator Newmark 0.5 0.25 \n');
    fprintf(fileID,'integrator HHT 0.9 \n');
    
    % Define analysis type
    fprintf(fileID,'analysis Transient \n');
elseif analysis.type == 2 || analysis.type == 3 % Pushover or Cyclic Analysis
    % Define Each Load Step
    [ control_node, control_dof, ~, step_size ] = fn_pushover_properties( first_story_node, analysis, story );
    int_controller = ['DisplacementControl ' num2str(control_node) ' ' num2str(control_dof) ' ' num2str(step_size)]; 
    fprintf(fileID,'integrator %s \n',int_controller);

    % Define analysis type
    fprintf(fileID,'analysis Static \n');
end

% Close File
fclose(fileID);

end

