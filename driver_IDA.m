%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Initial Setup
% Define Model
analysis.model_id = 6;
analysis.proceedure = 'NDP';
analysis.summit = 0;

% IDA Inputs
analysis.gm_seq_id = 10;
analysis.IDA_scale_factors = [0.5,0.7,0.9];
analysis.collapse_drift = 0.20;

% Secondary options
analysis.dead_load = 1;
analysis.live_load = 1;
analysis.opensees_SP = 1;
analysis.type = 1;
analysis.nonlinear = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.03;
analysis.hinge_stiff_mod = 10;
analysis.run_eigen = 0;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
% analysis.play_movie = 1;
% analysis.movie_scale = 1;

% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Import packages
import opensees.post_process.*
import opensees.write_tcl.*
import plotting_tools.fn_format_and_save_plot

% Sa Values
sa_x_t1 = 0.35;
sa_z_t1 = 0.88;

%% Define read and write directories
model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'model_data'];
tcl_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'opensees_data'];
outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'IDA'];
mkdir(outputs_dir)

% Load in Model Tables
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);

%% Run Opensees Model
for i = 1:length(analysis.IDA_scale_factors)
    collapse(i) = 0;
    
    % Change EQ Scale factor
    analysis.ground_motion_scale_factor = analysis.IDA_scale_factors(i);
    sa_x(i,1) = sa_x_t1*analysis.ground_motion_scale_factor;
    sa_z(i,1) = sa_z_t1*analysis.ground_motion_scale_factor;
    
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
    
    % Load ground motion data
    gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
    ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_seq_id,:);
    ground_motion_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
    if ground_motion_seq.eq_id_x ~= 0
        ground_motion.x = ground_motion_table(ground_motion_table.id == ground_motion_seq.eq_id_x,:);
    end
    if ground_motion_seq.eq_id_z ~= 0
        ground_motion.z = ground_motion_table(ground_motion_table.id == ground_motion_seq.eq_id_z,:);
    end
    if ground_motion_seq.eq_id_y ~= 0
        ground_motion.y = ground_motion_table(ground_motion_table.id == ground_motion_seq.eq_id_y,:);
    end
    
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
            collapse(i) = 1;
        else
            error('Opensees Failed')
        end
    end
    
    %% Post Process Data
    for n = 1:height(node)
        node_id = node.id(n);
           if node.record_disp(n)
               [ node_disp_raw ] = fn_xml_read([outputs_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
               node_disp_raw = node_disp_raw'; % flip to be node per row
               disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_x_TH']) = node_disp_raw(2,:);
               node.max_disp_x(n) = max(abs(node_disp_raw(2,:)));
               disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_z_TH']) = node_disp_raw(3,:);  % Currently Hard coded to three dimensions
               node.max_disp_z(n) = max(abs(node_disp_raw(3,:)));
           else
%                node.disp_x{n} = [];
               node.max_disp_x(n) = NaN;
%                node.disp_z{n} = [];
               node.max_disp_z(n) = NaN;
           end
    end
    
    % Calc Story Drift
    [ story.max_drift_x ] = fn_drift_profile( disp_TH, story, node, 'x' );
    [ story.max_drift_z ] = fn_drift_profile( disp_TH, story, node, 'z' );
    
    max_drift_x(i,1) = max(story.max_drift_x);
    max_drift_z(i,1) = max(story.max_drift_z);
end

%% Plot and Save Results
plot(max_drift_x,sa_x)
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( outputs_dir, 'IDA Plot EW Frame Direction', 2 )

plot(max_drift_z,sa_z)
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( outputs_dir, 'IDA Plot NS Wall Direction', 2 )

% Save results as csv
ida.sa_x = sa_x;
ida.sa_z = sa_z;
ida.max_drift_x = max_drift_x;
ida.max_drift_z = max_drift_z;
ida_table = struct2table(ida);
writetable(ida_table,[outputs_dir filesep 'ida_results_.csv'])
