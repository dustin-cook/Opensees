function [ ] = main_opensees_analysis( model, analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import opensees.*
import asce_41.*

% Create Write Directory
write_dir_opensees = ['outputs/' model.name{1} '/' analysis.proceedure '_' analysis.id '/opensees_data']; % TCL or Opensees does not like filesep command on windows, therefore must manually define forward slash seperators
if analysis.run_opensees % Don't clear the file if you don't want to run opensees
    fn_make_directory( write_dir_opensees )
end

% Define Read Directories
read_dir_model = [analysis.out_dir filesep 'model_data'];
read_dir_analysis = [analysis.out_dir filesep 'asce_41_data'];

% Load Model Data
node = readtable([read_dir_model filesep 'node.csv'],'ReadVariableNames',true);
element = readtable([read_dir_model filesep 'element.csv'],'ReadVariableNames',true);
story = readtable([read_dir_model filesep 'story.csv'],'ReadVariableNames',true);
joint = readtable([read_dir_model filesep 'joint.csv'],'ReadVariableNames',true);
hinge = readtable([read_dir_model filesep 'hinge.csv'],'ReadVariableNames',true);

% Define element hinge properties if not already defined
if analysis.nonlinear ~= 0 && analysis.model_type == 2 && ~exist([read_dir_analysis filesep 'element_analysis.mat'],'file')
    % Nonlinear MDOF with no capacities calculated
    ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
    [ element, joint ] = main_element_capacity( story, ele_props_table, element, analysis, joint, read_dir_analysis );
    [ element, ~ ] = main_hinge_properties( ele_props_table, element, joint );
    if ~exist(read_dir_analysis,'dir')
        fn_make_directory( read_dir_analysis )
    end
    save([read_dir_analysis filesep 'element_analysis.mat'],'element')
    save([read_dir_analysis filesep 'joint_analysis.mat'],'joint')
end

%% Run Opensees
% Define Number of Opensees runs to be performed
pushover_directions = {'x', '-x', 'z', '-z'};
if (analysis.type == 2 || analysis.type == 3) && strcmp(model.dimension,'3D') % 3D pushover or cyclic
    num_OS_runs = 4;
elseif (analysis.type == 2 || analysis.type == 3) && strcmp(model.dimension,'2D') % 2D pushover or cyclic
    num_OS_runs = 2;
else
    num_OS_runs = 1;
end

% Run for each specified run
for i = 1:num_OS_runs
    % Define inputs for this run
    analysis.pushover_direction = pushover_directions{i};

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
    
    % Define Hinge Group
    hinge_grouping = [];
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        hinge.group = zeros(height(hinge),1);
        hinge_grouping = 1:analysis.hinge_group_length:(height(hinge)+1);
        if hinge_grouping(end) < (height(hinge)+1)
            hinge_grouping = [hinge_grouping, height(hinge)+1];
        end
        for hg = 1:(length(hinge_grouping)-1)
            group_range = hinge.id >= hinge_grouping(hg) & hinge.id < hinge_grouping(hg+1);
            hinge.group(group_range) = hg;
        end
    end
    
    if analysis.run_opensees
        % Write TCL files
        main_write_tcl( model.dimension, write_dir_opensees, node, element, story, joint, hinge, analysis, read_dir_analysis, ground_motion, hinge_grouping );

        % Run Opensees
        fprintf('Running Opensees Analysis %i of %i \n',i,num_OS_runs)
        main_run_opensees( write_dir_opensees, analysis )
    end
end

% Postprocess OS data
main_post_process_opensees( analysis, model, story, node, element, joint, hinge, ground_motion, write_dir_opensees )

%% Write Summit Batch File
% if analysis.summit
%     write_dir_summit = [analysis.out_dir filesep 'summit'];
%     fn_make_directory( write_dir_summit )
%     fn_write_summit_batch_file( write_dir_summit, analysis.proceedure, model.name{1} )
% end
end

