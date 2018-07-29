%% INPUTS
model_name = 'ICBS_model_3D';
story_ht = [174 162 162 162 162 158];
story_dead_load = [2200 2325 2325 2325 2325 1900]*1000;
story_live_load = [128, 128, 128, 128, 128, 51]*1000;
bay_length.x = [71 300 300 300 300 300 71];
bay_length.z = [300 300 300];
damp_ratio = 0.05;
foundation = 'fix';
dims = '3D';
hazus_class_1 = 'C1';
hazus_class_2 = 'C2';

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
col_id =  [6  7  8  9  10 11 11 0  0  0  0  0  0];
beam_id = [12 12 12 12 13 14 15 16 17 0  0  0  0];
wall_id = [0  0  0  0  0  0  0  0  0  18 19 20 21];
x_offset_start = [0  0  0  0  0  0  0  0  0  12 60 0 0];
x_offset_end = [0  0  0  0  0  0  0  0  0  12 12 0 0];
grid_lines{1}.bays = 2:6;
grid_lines{2}.bays = 2:6;
grid_lines{3}.bays = 2:6;
grid_lines{4}.bays = 2:6;
grid_lines{5}.bays = 2:6;
grid_lines{6}.bays = 2:6;
grid_lines{7}.bays = 2:6;
grid_lines{8}.bays = 1;
grid_lines{9}.bays = 1;
grid_lines{10}.bays = 2;
grid_lines{11}.bays = 2;
grid_lines{12}.bays = 1:3;
grid_lines{13}.bays = 1:3;

direction = {'x' 'x' 'x' 'x' 'x' 'x' 'x' 'x' 'x' 'z' 'z' 'z' 'z'};

fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines, x_offset_start, x_offset_end)

%% Asseble Stories
story_group_id{1}.grid_id = [1 1 5 5 8 8 8 8 8 8 8 8 10 10 10 11];
story_group_id{1}.start_prim = [0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500 300 300 300 300];
story_group_id{1}.start_alt = [300 600 0 900 0 300 600 900 0 300 600 900 0 600 1200 900];
story_group_id{1}.direction = [1 1 1 1 1 1 1 1 1 1 1 1 3 3 3 3];

story_group_id{2}.grid_id = [2 2 6 6 9 9 9 9 9 9 9 9 12 12];
story_group_id{2}.start_prim = [0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500 0 0];
story_group_id{2}.start_alt = [300 600 0 900 0 300 600 900 0 300 600 900 -71 1571];
story_group_id{2}.direction = [1 1 1 1 1 1 1 1 1 1 1 1 3 3];

story_group_id{3}.grid_id = [3 3 6 6 9 9 9 9 9 9 9 9 13 13];
story_group_id{3}.start_prim = [0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500 0 0];
story_group_id{3}.start_alt = [300 600 0 900 0 300 600 900 0 300 600 900 -71 1571];
story_group_id{3}.direction = [1 1 1 1 1 1 1 1 1 1 1 1 3 3];

story_group_id{4}.grid_id = [4 4 6 6 9 9 9 9 9 9 9 9 13 13];
story_group_id{4}.start_prim = [0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500 0 0];
story_group_id{4}.start_alt = [300 600 0 900 0 300 600 900 0 300 600 900 -71 1571];
story_group_id{4}.direction = [1 1 1 1 1 1 1 1 1 1 1 1 3 3];

story_group_id{5}.grid_id = [4 4 7 7 9 9 9 9 9 9 9 9 13 13];
story_group_id{5}.start_prim = [0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500 0 0];
story_group_id{5}.start_alt = [300 600 0 900 0 300 600 900 0 300 600 900 -71 1571];
story_group_id{5}.direction = [1 1 1 1 1 1 1 1 1 1 1 1 3 3];

fn_assemble_story(input_dir, story_group_id);

%% Stack Stories
story_groups = [1 2 3 4 4 5];
fn_stack_stories(input_dir, story_ht, story_dead_load, story_live_load, story_groups, model_id)


