function [] = fn_main_IDA(analysis, model, story, node, gm_set_table, gm_idx, scale_factor, tcl_dir)
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
summary.collapse = 0;
analysis.ground_motion_scale_factor = scale_factor;
outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'IDA' '/' 'Scale_' num2str(analysis.ground_motion_scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];
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
fprintf(fileID,'puts "Defining Recorders Complete"\n');
fclose(fileID);

% Write Loads file
fn_define_loads( outputs_dir, analysis, node, model.dimension, story, 0, 0, ground_motion )

% Write Analysis Files
first_story_node = node.id(node.primary_story == 1);
fn_setup_analysis( outputs_dir, tcl_dir, analysis, first_story_node, story )
fn_define_analysis( outputs_dir, ground_motion, first_story_node, story.story_ht, analysis, story )

% Call Opensees
if analysis.summit
    command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' outputs_dir filesep 'run_analysis.tcl'];
else
    command = ['openseesSP ' outputs_dir filesep 'run_analysis.tcl'];
end
[status,cmdout] = system(command,'-echo');

% test for analysis failure and terminate Matlab
if status ~= 0
    if contains(cmdout,'Analysis Failure: Convergence') || contains(cmdout,'Analysis Failure: Singularity')
        warning('collapse or signularity')
        summary.collapse = 1;
    else
        error('Opensees Failed')
    end
end

clear cmdout

%% Post Process Data
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

% Calc Story Drift
[ story.max_drift_x ] = fn_drift_profile( disp_TH, story, node, 'x' );
[ story.max_drift_z ] = fn_drift_profile( disp_TH, story, node, 'z' );

clear disp_TH

summary.max_drift_x = max(story.max_drift_x);
summary.max_drift_z = max(story.max_drift_z);

% Save data for this run
save([outputs_dir filesep 'summary_results.mat'],'summary')
save([outputs_dir filesep 'story_analysis.mat'],'story')

end

