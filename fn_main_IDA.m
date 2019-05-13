function [ exist_status ] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gm_idx, scale_factor, tcl_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

import opensees.write_tcl.*
import opensees.post_process.*

% Defin gms for this run
ground_motion.x = gm_set_table(gm_idx,:);
ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};

% Iteration Parameters
analysis.ground_motion_scale_factor = scale_factor;
outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Scale_' num2str(analysis.ground_motion_scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];
mkdir(outputs_dir)

% Load spectral info and save Sa
model_periods = load([tcl_dir filesep 'model_analysis.mat']);
spectra_table = readtable([ground_motion.x.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
summary.sa_x = interp1(spectra_table.period,spectra_table.psa_5,model_periods.model.T1_x)*analysis.ground_motion_scale_factor;
clear spectra_table
spectra_table = readtable([ground_motion.z.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
summary.sa_z = interp1(spectra_table.period,spectra_table.psa_5,model_periods.model.T1_z)*analysis.ground_motion_scale_factor;
clear spectra_table

% Write Recorders File
file_name = [outputs_dir filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');
fprintf(fileID,'puts "Defining Recorders ..."\n');
fprintf(fileID,'setMaxOpenFiles 2000\n');
for n = 1:height(node)
   if node.record_disp(n)
        fprintf(fileID,'recorder Node -xml %s/nodal_disp_%s.xml -time -node %i -dof 1 3 disp\n',outputs_dir,num2str(node.id(n)),node.id(n));
   end
end
if analysis.nonlinear ~= 0 && ~isempty(hinge)
    for i = 1:height(hinge)
        hinge_y = node.y(node.id == hinge.node_1(i));
        if hinge_y == 0 && strcmp(hinge.direction{i},'primary')
            hinge_id = element.id(end) + hinge.id(i);
%             fprintf(fileID,'recorder Element %s %s/hinge_force_%s.%s -time -ele %s -dof 1 3 4 6 force \n', file_type, write_dir, num2str(hinge_id), file_ext, num2str(hinge_id));
            fprintf(fileID,'recorder Element -xml %s/hinge_deformation_%s.xml -time -ele %s deformation \n', outputs_dir, num2str(hinge_id), num2str(hinge_id));
        end
    end
end
fclose(fileID);

% Write Loads file
fn_define_loads( outputs_dir, analysis, node, model.dimension, story, 0, 0, ground_motion )

% Write Analysis Files
first_story_node = node.id(node.primary_story == 1);
fn_setup_analysis( outputs_dir, tcl_dir, analysis, first_story_node, story )
fn_define_analysis( outputs_dir, ground_motion, first_story_node, story.story_ht, analysis, story )

% Call Opensees
fprintf('Running Opensess... \n')
if analysis.summit
    command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' outputs_dir filesep 'run_analysis.tcl'];
else
    command = ['openseesSP ' outputs_dir filesep 'run_analysis.tcl'];
end
if analysis.suppress_outputs
    [status,cmdout] = system(command);
else
    [status,cmdout] = system(command,'-echo');
end

% test for analysis failure and terminate Matlab
exist_status = 0;
if contains(cmdout,'Analysis Failure: Collapse')
    summary.collapse = 1; % Collapse triggered by drift limit
    fprintf('Model Reached Collapse Limit \n')
elseif contains(cmdout,'Analysis Failure: Convergence') || contains(cmdout,'Analysis Failure: Singularity')
    summary.collapse = 2; % Collapse triggered by convergence or singularity issue
    fprintf('Unexpected Opensees failure \n')
    fprintf('Model Experienced a Convergence Failure (Treat as collapsed)')
elseif status ~= 0
    summary.collapse = 3; % Unexpected Opensees failure (shouldnt get here)
    fprintf('UNHANDLED OPENSEES FAILURE \n')
    exist_status = 1;
else
    summary.collapse = 0;
    fprintf('Model Ran Successfully \n')
end

fprintf('Opensees Completed \n')
clear cmdout

%% Post Process Data
if exist_status
    fprintf('Skipping Postprocessing of Opensees')
else
    fprintf('Postprocessing Opensees Ouputs from Directory: %s \n',outputs_dir)
    % Nodal displacements
    for n = 1:height(node)
        node_id = node.id(n);
           if node.record_disp(n)
               [ node_disp_raw ] = fn_xml_read([outputs_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
               node_disp_raw = node_disp_raw'; % flip to be node per row
               disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_x_TH']) = node_disp_raw(2,:);
               disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_z_TH']) = node_disp_raw(3,:);  % Currently Hard coded to three dimensions
           end
           clear node_disp_raw
    end

    % Hinge Deformations
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        for i = 1:height(hinge)
            hinge_y = node.y(node.id == hinge.node_1(i));
            if hinge_y == 0 && strcmp(hinge.direction{i},'primary')
                hinge_id = element.id(end) + hinge.id(i);
                [ hinge_deformation_TH ] = fn_xml_read([outputs_dir filesep 'hinge_deformation_' num2str(hinge_id) '.xml']);
                hinge.max_deform(i) = max(abs(hinge_deformation_TH(:,2)));
                clear hinge_deformation_TH
            else
                hinge.max_deform(i) = NaN;
            end
        end
    end
    save([outputs_dir filesep 'hinge_analysis.mat'],'hinge')

    % Calc Story Drift
    [ story.max_drift_x ] = fn_drift_profile( disp_TH, story, node, 'x' );
    [ story.max_drift_z ] = fn_drift_profile( disp_TH, story, node, 'z' );

    clear disp_TH

    summary.max_drift_x = max(story.max_drift_x);
    summary.max_drift_z = max(story.max_drift_z);

    % Save data for this run
    fprintf('Writing IDA Results to Directory: %s \n',outputs_dir)
    save([outputs_dir filesep 'summary_results.mat'],'summary')
    save([outputs_dir filesep 'story_analysis.mat'],'story')
end

end

