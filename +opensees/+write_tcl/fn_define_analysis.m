function [ ] = fn_define_analysis( output_dir, ground_motion, first_story_node, story_ht, analysis )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Write Loads File
file_name = [output_dir filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

%% Run the Analysis
if analysis.type == 1 % Dynamic
    % Set analysis timestep equal to the ground motion time step
    time_step = ground_motion.x.eq_dt; % Clean up to be not based on x?
    num_steps = ground_motion.x.eq_length;

    fprintf(fileID,'source %s/setup_dynamic_analysis.tcl \n', output_dir);
    fprintf(fileID,'set singularity_check 0 \n');
    fprintf(fileID,'set collapse_check 0 \n');
    fprintf(fileID,'set currentTime [getTime] \n');
    fprintf(fileID,'set ok 0 \n');

    if analysis.solution_algorithm
        % While loop through each step of the ground motion
        fprintf(fileID,'while {$ok == 0 && $currentTime < %f && $collapse_check == 0 && $singularity_check == 0} { \n',time_step*num_steps);
        tolerance = 1e-6;
        fprintf(fileID,'test NormDispIncr %f 100 \n',tolerance);
        fprintf(fileID,'algorithm KrylovNewton \n');
        fprintf(fileID,'set dt_reduce %f \n', 1);
        fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
        fprintf(fileID,'set dt_min [expr %f/($dt_reduce*100)] \n', time_step);
        fprintf(fileID,'set ok [analyze 1 $dt $dt_min $dt] \n');
        % fprintf(fileID,'set ok [analyze %i $dt] \n',round(num_steps));
        fprintf(fileID,'puts "analysis failure = $ok " \n');

        % Loop Through Algorithms
        algorithm_typs = {'NewtonLineSearch', 'Newton -initial', 'Newton'};
        for a = 1:length(algorithm_typs)
            fprintf(fileID,'if {$ok != 0} { \n');
            fprintf(fileID,'puts "analysis failed, try %s" \n', algorithm_typs{a});
        %     fprintf(fileID,'source %s/setup_dynamic_analysis.tcl \n', output_dir);
            fprintf(fileID,'test NormDispIncr %f 100 \n',tolerance);
            fprintf(fileID,'algorithm %s \n', algorithm_typs{a});
            fprintf(fileID,'set dt_reduce %f \n', 1);
            fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
            fprintf(fileID,'set ok [analyze 1 $dt] \n');
            fprintf(fileID,'} \n');
        end

        % Loop Though dt
        for t = 1:8
            fprintf(fileID,'if {$ok != 0} { \n');
            dt_reduction = 10*t;
            fprintf(fileID,'puts "analysis failed, try dt/%f" \n', dt_reduction);
        %     fprintf(fileID,'source %s/setup_dynamic_analysis.tcl \n', output_dir);
            fprintf(fileID,'test NormDispIncr %f 100 \n',tolerance);
            fprintf(fileID,'set dt_reduce %f \n', dt_reduction);
            fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
            fprintf(fileID,'set ok [analyze 1 $dt] \n');
            fprintf(fileID,'} \n');
        end

        % Loop Through Tolerance
        for tol = 1:4
            fprintf(fileID,'if {$ok != 0} { \n');
            tolerance = 1e-5*10^(tol);
            fprintf(fileID,'puts "analysis failed, try tolerance = %f" \n', tolerance);
        %     fprintf(fileID,'source %s/setup_dynamic_analysis.tcl \n', output_dir);
            fprintf(fileID,'test NormDispIncr %f 1000 \n', tolerance);
            fprintf(fileID,'set dt_reduce %f \n', dt_reduction);
            fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
            fprintf(fileID,'set ok [analyze 1 $dt] \n');
            fprintf(fileID,'} \n');
        end

        fprintf(fileID,'set currentTime [getTime] \n');
        fprintf(fileID,'puts "Time = $currentTime" \n');

        % Check for singularity and collapse
        fprintf(fileID,'if {$ok == 0} { \n');
        % Define Displacement
        fprintf(fileID,'set node_at_floor_1 %i \n', first_story_node(1));
        fprintf(fileID,'set floor_displ_1 "[nodeDisp $node_at_floor_1 1]" \n');
        fprintf(fileID,'puts "First Story Disp = $floor_displ_1" \n');
        fprintf(fileID,'set height_floor_1 %f \n', story_ht(1));
        fprintf(fileID,'set floor_drift_1 [expr abs($floor_displ_1/$height_floor_1)] \n');
        % Check for Singularity
        fprintf(fileID,'set check_QNAN_1 [string first QNAN $floor_displ_1 1] \n');
        fprintf(fileID,'puts "QNAN Check = $check_QNAN_1" \n');
        fprintf(fileID,'set check_IND_1 [string first IND $floor_displ_1 1] \n');
        fprintf(fileID,'puts "IND Check = $check_IND_1" \n');
        fprintf(fileID,'if {($floor_displ_1 > 1000000) || ($check_QNAN_1 != -1) || ($check_IND_1 != -1)} { \n');
        fprintf(fileID,'set singularity_check 1 \n');
        fprintf(fileID,'} \n');
        fprintf(fileID,'puts "Singularity = $singularity_check" \n');
        % Check for Collapse
        if analysis.collapse_drift > 0
            fprintf(fileID,'if {$floor_drift_1 > %f} { \n', analysis.collapse_drift);
            fprintf(fileID,'set collapse_check 1 \n');
            fprintf(fileID,'} \n');
            fprintf(fileID,'puts "Collapse = $collapse_check" \n');
        end
        fprintf(fileID,'} \n');
        fprintf(fileID,'} \n');
    else
        tolerance = 1e-6;
        fprintf(fileID,'test NormDispIncr %f 100 \n',tolerance);
        fprintf(fileID,'algorithm Newton \n');
        fprintf(fileID,'set dt_reduce %f \n', analysis.initial_timestep_factor);
        fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
        fprintf(fileID,'puts "Analysis dt = $dt" \n');
        fprintf(fileID,'set ok [analyze %i $dt] \n',round(num_steps*analysis.initial_timestep_factor));
        fprintf(fileID,'puts "analysis failure = $ok " \n');
    end
elseif analysis.type == 2 % Pushover
    fprintf(fileID,'source %s/setup_pushover_analysis.tcl \n', output_dir);
    fprintf(fileID,'set singularity_check 0 \n');
    fprintf(fileID,'set collapse_check 0 \n');
    fprintf(fileID,'set currentTime [getTime] \n');
    fprintf(fileID,'set ok 0 \n');
    fprintf(fileID,'set ok [analyze %i] \n', analysis.pushover_num_steps);
    fprintf(fileID,'puts "analysis failure = $ok " \n');
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
    fprintf(fileID,'set time_step_file %s/final_time_step_reduction.txt \n',output_dir);
    fprintf(fileID,'set TS [open $time_step_file "w"] \n');
    fprintf(fileID,'puts $TS " $dt_reduce" \n');
    fprintf(fileID,'close $TS \n');
end
fprintf(fileID,'wipe \n');
fprintf(fileID,'} \n');

% Close File
fclose(fileID);


end

