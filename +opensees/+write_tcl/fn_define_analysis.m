function [ ] = fn_define_analysis( output_dir, ground_motion, first_story_node, story_ht, analysis, story )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Write Loads File
file_name = [output_dir filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

%% Set Parameters
min_tolerance_steps = 10;
max_tolerance_steps = 1000;

%% Initial Analysis Setup
fprintf(fileID,'source %s/setup_analysis.tcl \n', output_dir);
fprintf(fileID,'set singularity_check 0 \n');
fprintf(fileID,'set collapse_check 0 \n');
fprintf(fileID,'set currentTime [getTime] \n');
fprintf(fileID,'set ok 0 \n');

%% Run the Analysis
if analysis.type == 1 % Dynamic 
    % Set analysis timestep equal to the ground motion time step
    time_step = ground_motion.x.eq_dt; % Clean up to be not based on x?
    num_steps = ground_motion.x.eq_length;
    
    % Set loop factors
    dt_reduction = [1,10,20,50,100];
    algorithm_typs = {'KrylovNewton', 'NewtonLineSearch', 'Newton -initial', 'Newton'};
    tolerance = [1e-5, 1e-4 0.001, 0.01, 0.1, 1];

    if analysis.solution_algorithm
        % Set up Log Files
        log_file = [output_dir '/converge_tol_file.txt'];
        fprintf(fileID,'set converge_tol_file [open %s w] \n', log_file);
        convergence_file = [output_dir '/converge_file.txt'];
        fprintf(fileID,'set converge_file [open %s w] \n', convergence_file);
    
        % While loop through each step of the ground motion
        eq_total_time = time_step*num_steps;
        fprintf(fileID,'while {$ok == 0 && $currentTime < %f && $collapse_check == 0 && $singularity_check == 0} { \n',eq_total_time);
        
        % Output current time progress
        fprintf(fileID,'puts "tFinal is %f; and tCurrent is $currentTime" \n',eq_total_time);
        
        % Run analysis with basic props
        fprintf(fileID,'set tol %f \n', tolerance(1));
        fprintf(fileID,'test NormDispIncr $tol %i \n', min_tolerance_steps);
        fprintf(fileID,'algorithm KrylovNewton \n');
        fprintf(fileID,'set dt_reduce %f \n', 1);
        fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
        fprintf(fileID,'set dt_max %f \n', time_step);
        fprintf(fileID,'set dt_min [expr %f/($dt_reduce*100)] \n', time_step);
        fprintf(fileID,'set ok [analyze 1 $dt $dt_min $dt_max] \n');
        % fprintf(fileID,'set ok [analyze %i $dt] \n',round(num_steps));
        fprintf(fileID,'puts "analysis failure = $ok " \n');

        % Loop Through Tolerance
        for tol = 1:length(tolerance)
            
            % Loop Though dt
            for t = 1:length(dt_reduction)

                % Loop Through Algorithms
%                 for a = 1:length(algorithm_typs)
                    fprintf(fileID,'if {$ok != 0} { \n');
                    fprintf(fileID,'puts "analysis failed, try try tolerance = %f, dt/%f, and %s" \n', tolerance(tol), dt_reduction(t), algorithm_typs{1});
                    fprintf(fileID,'set tol %f \n', tolerance(tol));
                    if tol <= 4
                        fprintf(fileID,'test NormDispIncr $tol %i \n', min_tolerance_steps);
                    else
                        fprintf(fileID,'test NormDispIncr $tol %i \n', max_tolerance_steps);
                    end
                    fprintf(fileID,'algorithm %s \n', algorithm_typs{1});
                    fprintf(fileID,'set dt_reduce %f \n', dt_reduction(t));
                    fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
                    fprintf(fileID,'set ok [analyze 1 $dt] \n');
                    fprintf(fileID,'} \n');
%                 end
            end
        end

        % Output time
        fprintf(fileID,'set currentTime [getTime] \n');
        fprintf(fileID,'puts "Time = $currentTime" \n');
        
        % Save info to log file
        fprintf(fileID,'set converge_tol_log "$currentTime $tol" \n');
        fprintf(fileID,'puts $converge_tol_file $converge_tol_log \n');

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
        
        % Check if Model Fully Converged
        fprintf(fileID,'if {$currentTime > %f || $singularity_check == 1 || $collapse_check == 1} { \n',eq_total_time-2); % 2 seconds before the end of the EQ is okay
        fprintf(fileID,'puts "EQ fully converged" \n');
        fprintf(fileID,'puts $converge_file 1 \n');
        fprintf(fileID,'} else { \n');
        fprintf(fileID,'puts "EQ NOT FULLY converged" \n');
        fprintf(fileID,'puts $converge_file 0 \n');
        fprintf(fileID,'} \n');
    
        % Close Log Files
        fprintf(fileID,'close $converge_tol_file \n');
        fprintf(fileID,'close $converge_file \n');
    else
        fprintf(fileID,'set dt_reduce %f \n', analysis.initial_timestep_factor);
        fprintf(fileID,'set dt [expr %f/$dt_reduce] \n', time_step);
        fprintf(fileID,'puts "Analysis dt = $dt" \n');
        fprintf(fileID,'set ok [analyze %i $dt] \n',round(num_steps*analysis.initial_timestep_factor));
        fprintf(fileID,'puts "analysis failure = $ok " \n');
    end
elseif analysis.type == 2 % Pushover
    fprintf(fileID,'set ok [analyze %i] \n', analysis.pushover_num_steps);
    fprintf(fileID,'puts "analysis failure = $ok " \n');
elseif analysis.type == 3 % Static Cyclic
    fprintf(fileID,'source %s/setup_static_cyclic_analysis.tcl \n', output_dir);
    fprintf(fileID,'set singularity_check 0 \n');
    fprintf(fileID,'set collapse_check 0 \n');
    fprintf(fileID,'set currentTime [getTime] \n');
    fprintf(fileID,'set ok 0 \n');
    % Define Each Load Step
    control_node = first_story_node(end);
    if strcmp(analysis.pushover_direction,'x')
        control_dof = 1;
    elseif strcmp(analysis.pushover_direction,'z')
        control_dof = 3;
    end
    max_displacement = analysis.pushover_drift*(story.y_start(end)+story.story_ht(end));
    step_size = max_displacement / analysis.pushover_num_steps;
    steps = ones(1,analysis.pushover_num_steps)*step_size;
    cycle_cuts = round((ones(1,4)*length(steps)).*[1/4, 1/2, 3/4, 1]);
    load_pattern = [];
    for i = 1:4
        cycle =  [steps(1:cycle_cuts(i)),-steps(1:cycle_cuts(i)),-steps(1:cycle_cuts(i)),steps(1:cycle_cuts(i))];
        load_pattern = [load_pattern, cycle, cycle, cycle];
    end
    % Run through each step of load pattern
    for i = 1:length(load_pattern)
        fprintf(fileID,'integrator DisplacementControl %i %i %f \n',control_node,control_dof,load_pattern(i));
        fprintf(fileID,'analysis Static \n');
        fprintf(fileID,'set ok [analyze 1] \n');
    end
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

