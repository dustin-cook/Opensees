clear
close
clc

%% INPUTS
model_name = 'ICBS_model_5ew';
story_ht = [174 162 162 162 162 158];
story_dead_load = (1/2)*[2200 2325 2325 2325 2325 1900]*1000;
story_live_load = (1/2)*[128, 128, 128, 128, 128, 51]*1000;
bay_length.x = [300 300 300 300 300 100 300 300 300 300 300];
bay_length.z = [0];
damp_ratio = 0.05;
foundation = 'fix';
dims = '2D';
hazus_class_1 = 'C1';
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
    model_table.dimension = dims;
    model_table.hazus_class_1 = hazus_class_1;
    model_table.hazus_class_2 = hazus_class_2;
    writetable(model_table,['inputs' filesep 'model.csv'],'WriteVariableNames',true);
else
    model_id = model_table.id(checker);
end

%% Build Grid Lines
col_id = [6 7 8 9 0 10 11 11];
beam_id = [12 12 12 12 1 13 14 15];
wall_id = [0 0 0 0 0 0 0 0];
grid_lines{1}.bays = 7:11;
grid_lines{2}.bays = 7:11;
grid_lines{3}.bays = 7:11;
grid_lines{4}.bays = 7:11;
grid_lines{5}.bays = 6;
grid_lines{6}.bays = 1:5;
grid_lines{7}.bays = 1:5;
grid_lines{8}.bays = 1:5;

direction = {'x' 'x' 'x' 'x' 'x' 'x' 'x' 'x'};

fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines)

%% Asseble Stories
story_group_id{1}.grid_id = [6 1 5];
story_group_id{1}.start_prim = [0 1600 1500];
story_group_id{1}.start_alt = [0 0 0];
story_group_id{1}.direction = [1 1 1];

story_group_id{2}.grid_id = [7 2 5];
story_group_id{2}.start_prim = [0 1600 1500];
story_group_id{2}.start_alt = [0 0 0];
story_group_id{2}.direction = [1 1 1];

story_group_id{3}.grid_id = [7 3 5];
story_group_id{3}.start_prim = [0 1600 1500];
story_group_id{3}.start_alt = [0 0 0];
story_group_id{3}.direction = [1 1 1];

story_group_id{4}.grid_id = [7 4 5];
story_group_id{4}.start_prim = [0 1600 1500];
story_group_id{4}.start_alt = [0 0 0];
story_group_id{4}.direction = [1 1 1];

story_group_id{5}.grid_id = [8 4 5];
story_group_id{5}.start_prim = [0 1600 1500];
story_group_id{5}.start_alt = [0 0 0];
story_group_id{5}.direction = [1 1 1];


fn_assemble_story(input_dir, story_group_id);

%% Stack Stories
story_groups = [1 2 3 4 4 5];
fn_stack_stories(input_dir, story_ht, story_dead_load, story_live_load, story_groups, model_id)


