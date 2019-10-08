function [ exit_status ] = fn_main_IDA(analysis, model, story, element, node, hinge, joint, ground_motion, scale_factor, building_period, tcl_dir, ida_opensees_dir, ida_summary_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

import opensees.write_tcl.*

% Load spectral info and save Sa
spectra_table = readtable([ground_motion.x.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
summary.sa_x = interp1(spectra_table.period,spectra_table.psa_5,building_period(1))*scale_factor;
if analysis.run_z_motion
    spectra_table = readtable([ground_motion.z.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
    summary.sa_z = interp1(spectra_table.period,spectra_table.psa_5,building_period(2))*scale_factor;
end

% Write Recorders File
fn_define_recorders( ida_opensees_dir, model.dimension, node, element, joint, hinge, analysis )

% Write Loads file
analysis.ground_motion_scale_factor = scale_factor;
fn_define_loads( ida_opensees_dir, analysis, node, model.dimension, story, element, ground_motion )

% Write Analysis Files
first_story_node = node.id(node.primary_story == 1);
fn_setup_analysis( ida_opensees_dir, tcl_dir, analysis, first_story_node, story )
fn_define_analysis( ida_opensees_dir, ground_motion, first_story_node, story.story_ht, analysis, story )

% Call Opensees
fprintf('Running Opensess... \n')
if analysis.summit
    command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' ida_opensees_dir filesep 'run_analysis.tcl'];
else
    command = ['openseesSP ' ida_opensees_dir filesep 'run_analysis.tcl'];
end
if analysis.suppress_outputs
    [status,cmdout] = system(command);
else
    [status,cmdout] = system(command,'-echo');
end

% test for analysis failure and terminate Matlab
exit_status = 0;
if contains(cmdout,'Analysis Failure: Collapse')
    summary.collapse = 1; % Collapse triggered by drift limit
    fprintf('Model Reached Collapse Limit \n')
elseif contains(cmdout,'Analysis Failure: Singularity')
    summary.collapse = 2; % Collapse triggered by singularity issue
    fprintf('Unexpected Opensees failure \n')
    fprintf('Model Experienced a Singularity Failure (Treat as collapsed)')
elseif contains(cmdout,'Analysis Failure: Convergence')
    summary.collapse = 3; % Collapse triggered by convergence
    fprintf('Unexpected Opensees failure \n')
    fprintf('Model Experienced a Convergence Failure (Treat as collapsed)')
elseif status ~= 0
    summary.collapse = 4; % Unexpected Opensees failure (shouldnt get here)
    fprintf('UNHANDLED OPENSEES FAILURE \n')
    exit_status = 1;
else
    summary.collapse = 0;
    fprintf('Model Ran Successfully \n')
end

fprintf('Opensees Completed \n')

% Save summary data
save([ida_summary_dir filesep 'summary_results.mat'],'summary')
end

