function [ ] = main_opensees_analysis( model, analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import opensees.*

% Create Write Directory
write_dir_opensees = ['outputs/' model.name{1} '/' analysis.proceedure '/opensees_data']; % TCL or Opensees does not like filesep command on windows, therefore must manually define forward slash seperators
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

%% Run Opensees
if analysis.run_opensees
    % Define Number of Opensees runs to be performed
    pushover_directions = {'x', 'z'};
    if analysis.type == 2 && strcmp(model.dimension,'3D') % 3D pushover
        num_OS_runs = 2;
    else
        num_OS_runs = 1;
    end

    for i = 1:num_OS_runs
        % Define inputs for this run
        analysis.pushover_direction = pushover_directions{i};

        % Write TCL files
        [ node, ground_motion ] = main_write_tcl( model.dimension, write_dir_opensees, node, element, story, joint, hinge, analysis, read_dir_analysis );

        % Run Opensees
        main_run_opensees( write_dir_opensees )
    end
    % Postprocess OS data
    main_post_process_opensees( analysis, model, story, node, element, hinge, ground_motion, write_dir_opensees )
end

%% Write Summit Batch File
if analysis.summit_SP
    write_dir_summit = [analysis.out_dir filesep 'summit'];
    fn_make_directory( write_dir_summit )
    fn_write_summit_batch_file( write_dir_summit, analysis.proceedure, model.name{1} )
end
end

