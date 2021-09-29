function [ exit_status ] = fn_main_IDA(analysis, model, story, element, node, hinge, joint, ground_motion, tcl_dir, ida_opensees_dir, ida_summary_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

import opensees.write_tcl.*

% Load spectral info and save Sa
summary.sa_x = ground_motion.x.sa;
summary.pga_x = analysis.ground_motion_scale_factor * ground_motion.x.pga;
if analysis.run_z_motion
    summary.sa_z = ground_motion.z.sa;
    summary.pga_z = analysis.ground_motion_scale_factor * ground_motion.z.pga;
end

% Write Recorders File
if analysis.general_ida
    fn_define_recorders_ida( ida_opensees_dir, model.dimension, node.id, element )
else
    fn_define_recorders( ida_opensees_dir, model.dimension, node, element, joint, hinge, analysis )
end


% Write Loads file
if analysis.general_ida
    fn_define_loads_ida( ida_opensees_dir, analysis.ground_motion_scale_factor, model.dimension, ground_motion, analysis.g_unit )
else
    fn_define_loads( ida_opensees_dir, analysis, node, model.dimension, story, element, joint, ground_motion, model )
end

% Write Analysis Files
first_story_node = node.id(node.primary_story == 1);
if analysis.general_ida
    fn_setup_analysis( ida_opensees_dir, tcl_dir, analysis )
    fn_define_analysis( ida_opensees_dir, ground_motion, first_story_node, story.story_ht, analysis )
else
    fn_setup_analysis( ida_opensees_dir, tcl_dir, analysis, first_story_node, story )
    fn_define_analysis( ida_opensees_dir, ground_motion, first_story_node, story.story_ht, analysis, story )
end

% Call Opensees
fprintf('Running Opensess... \n')
if analysis.summit
    command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' ida_opensees_dir filesep 'run_analysis.tcl'];
else
    if analysis.opensees_SP
        command = ['openseesSP ' ida_opensees_dir filesep 'run_analysis.tcl'];
    else
        command = ['opensees ' ida_opensees_dir filesep 'run_analysis.tcl'];
    end
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
    summary.collapse = 0; % Didn't catch anything (also shouldnt get here)
    fprintf('Model Ran Successfully \n')
end

fprintf('Opensees Completed \n')

%% Save summary data
save([ida_summary_dir filesep 'summary_results.mat'],'summary')

% % Copy primary node files to the summary driver
% for i = 1:length(node.primary_nodes)
%     try
%         source = [ida_opensees_dir filesep 'nodal_disp_' num2str(node.primary_nodes(i)) '.*'];
%         destination = ida_summary_dir;
%         copyfile(source,destination)
%     catch
%         warning('No matching files were found.')
%     end
% end

end

