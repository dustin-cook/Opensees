clear
close
clc

%% INPUTS
model_name = 'rob_sdof';
story_ht = [240];
story_dead_load = [48]*1000;
story_live_load = [0]*1000;
bay_length.x = [0];
bay_length.z = [0];
damp_ratio = 0.05;
foundation = 'fix';
dims = '2D';
hazus_class_1 = 'NA';
hazus_class_2 = 'NA';

%% Import Packages
import tools.*

%% Create Input Directory
input_dir = ['inputs' filesep 'models' filesep model_name];
if ~exist(input_dir,'dir')
    mkdir(input_dir);
end

%% Load and Write to Model Database
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
checker = strcmp(model_table.name,model_name);
if sum(checker) == 0
    model_id = length(model_table.id) + 1;
    model_table.id(model_id) = model_id;
    model_table.name{model_id} = model_name;
    model_table.foundation{model_id} = foundation;
    model_table.damp_ratio(model_id) = damp_ratio;
    model_table.dimension{model_id} = dims;
    model_table.hazus_class_1{model_id} = hazus_class_1;
    model_table.hazus_class_2{model_id} = hazus_class_2;
    writetable(model_table,['inputs' filesep 'model.csv'],'WriteVariableNames',true);
else
    model_id = model_table.id(checker);
end

%% Build Grid Lines
col_id = [23];
beam_id = [0];
wall_id = [0];
x_offset_start = [0];
x_offset_end = [0];
trib_wt = [1];

grid_lines{1}.bays = 1;

direction = {'x'};

fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines, x_offset_start, x_offset_end, trib_wt)

%% Asseble Stories
story_group_id{1}.grid_id = [1];
story_group_id{1}.start_prim = [0];
story_group_id{1}.start_alt = [0];
story_group_id{1}.direction = [1];

fn_assemble_story(input_dir, story_group_id);

%% Stack Stories
story_groups = [1];
fn_stack_stories(input_dir, story_ht, story_dead_load, story_live_load, story_groups, model_id)


