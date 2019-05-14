function [ ] = fn_define_analysis( write_dir, ground_motion, primary_nodes, story_ht, analysis, story )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import opensees.write_tcl.*

%% Write Loads File
file_name = [write_dir filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

%% Initial Analysis Setup
fprintf(fileID,'source %s/setup_analysis.tcl \n', write_dir);
fprintf(fileID,'set singularity_check 0 \n');
fprintf(fileID,'set collapse_check 0 \n');
fprintf(fileID,'set currentStep [getTime] \n');
fprintf(fileID,'set ok 0 \n');

%% Run the Analysis
fprintf(fileID,'puts "Analyzing Model ..." \n');

if analysis.type == 1 % Dynamic 
    % Set analysis timestep equal to the ground motion time step
    time_step = ground_motion.x.eq_dt; % Clean up to be not based on x?
    num_steps = ground_motion.x.eq_length;

    if analysis.solution_algorithm
        % Solution Algorithm Setup
        eq_total_time = time_step*num_steps;
        fn_solution_algorithm( fileID, analysis, write_dir, eq_total_time, time_step, primary_nodes, story_ht )
    else
        fprintf(fileID,'set dt_reduce %f \n', analysis.initial_timestep_factor);
        fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
        fprintf(fileID,'puts "Analysis dt = $dt" \n');
        fprintf(fileID,'set ok [analyze %i $dt] \n',round(num_steps*analysis.initial_timestep_factor));
        fprintf(fileID,'puts "analysis failure = $ok " \n');
    end
elseif analysis.type == 2 || analysis.type == 3 % Pushover or Static Cyclic
    if analysis.solution_algorithm
        
        % Calc other pushover analysis properties
        [ control_node, control_dof, max_displacement, step_size ] = fn_pushover_properties( primary_nodes, analysis, story );
        
        % Call solution algorithm
        fn_solution_algorithm( fileID, analysis, write_dir, max_displacement, step_size, primary_nodes, story_ht, control_node, control_dof )
    else
        fprintf(fileID,'set ok [analyze %i] \n', analysis.pushover_num_steps);
        fprintf(fileID,'puts "analysis failure = $ok " \n');
    end
end

%% Specify Analysis Failure or Success
% Convergence Failure
fprintf(fileID,'if {$ok != 0} { \n');
fprintf(fileID,'puts "Analysis Failure: Convergence" \n');
fprintf(fileID,'wipe \n');
fprintf(fileID,'exit \n');
fprintf(fileID,'} \n');
% Singularity Failure
fprintf(fileID,'if {$singularity_check == 1} { \n');
fprintf(fileID,'puts "Analysis Failure: Singularity" \n');
fprintf(fileID,'wipe \n');
fprintf(fileID,'exit \n');
fprintf(fileID,'} \n');
% Collapse Failure
fprintf(fileID,'if {$collapse_check == 1} { \n');
fprintf(fileID,'puts "Analysis Failure: Collapse" \n');
fprintf(fileID,'wipe \n');
fprintf(fileID,'exit \n');
fprintf(fileID,'} \n');
% Analysis Success
fprintf(fileID,'if {$ok == 0} { \n');
fprintf(fileID,'puts "Analysis Complete!" \n');
% Write output file for final analysis time step
if analysis.type == 1
    fprintf(fileID,'set time_step_file %s/final_time_step_reduction.txt \n',write_dir);
    fprintf(fileID,'set TS [open $time_step_file "w"] \n');
    fprintf(fileID,'puts $TS " $dt_reduce" \n');
    fprintf(fileID,'close $TS \n');
end
fprintf(fileID,'wipe \n');
fprintf(fileID,'exit \n');
fprintf(fileID,'} \n');

% Close File
fclose(fileID);


end

