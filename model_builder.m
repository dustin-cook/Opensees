%% INPUTS
model_name = 'ICBS_model_3D';
story_ht = [174 162 162 162 162 158];
story_wt = [2200 2325 2325 2325 2325 1900];
bay_length.x = [71 300 300 300 300 300 71];
bay_length.z = [300 300 300];
damp_ratio = 0.05;
foundation = 'fix';

%% Import Packages
import tools.*

%% Create Input Directory
input_dir = ['inputs/' model_name];
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
    writetable(model_table,['inputs' filesep 'model.csv'],'WriteVariableNames',true);
else
    model_id = model_table.id(checker);
end

%% Build Grid Lines
col_id = [6 0 6 7 0 0 0 7 0];
beam_id = [8 0 9 10 0 16 17 11 0];
wall_id = [0 12 0 0 14 0 0 0 15];
grid_lines{1}.bays = 2:6;
grid_lines{2}.bays = 1;
grid_lines{3}.bays = 2:6;
grid_lines{4}.bays = 2:6;
grid_lines{5}.bays = 1:3;
grid_lines{6}.bays = 1;
grid_lines{7}.bays = 1;
grid_lines{8}.bays = 2:6;
grid_lines{9}.bays = 1:3;

direction = {'x' 'z' 'x' 'x' 'z' 'x' 'x' 'x' 'z'};

fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines)

%% Asseble Stories
story_group_id{1}.grid_id = [1 1 1 1 2 2 2 2 6 6 6 6 6 6 6 6];
story_group_id{1}.start_prim = [0 0 0 0 300 300 300 300 -71 -71 -71 -71 1500 1500 1500 1500];
story_group_id{1}.start_alt = [0 300 600 900 0 600 900 1200 0 300 600 900 0 300 600 900];
story_group_id{1}.direction = [1 1 1 1 3 3 3 3 1 1 1 1 1 1 1 1];

story_group_id{2}.grid_id = [4 3 3 4 5 5 7 7 7 7 7 7 7 7];
story_group_id{2}.start_prim = [0 0 0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500];
story_group_id{2}.start_alt = [0 300 600 900 -71 1571 0 300 600 900 0 300 600 900];
story_group_id{2}.direction = [1 1 1 1 3 3 1 1 1 1 1 1 1 1];

story_group_id{3}.grid_id = [4 3 3 4 9 9 7 7 7 7 7 7 7 7];
story_group_id{3}.start_prim = [0 0 0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500];
story_group_id{3}.start_alt = [0 300 600 900 -71 1571 0 300 600 900 0 300 600 900];
story_group_id{3}.direction = [1 1 1 1 3 3 1 1 1 1 1 1 1 1];

story_group_id{4}.grid_id = [8 3 3 8 9 9 7 7 7 7 7 7 7 7];
story_group_id{4}.start_prim = [0 0 0 0 0 0 -71 -71 -71 -71 1500 1500 1500 1500];
story_group_id{4}.start_alt = [0 300 600 900 -71 1571 0 300 600 900 0 300 600 900];
story_group_id{4}.direction = [1 1 1 1 3 3 1 1 1 1 1 1 1 1];

fn_assemble_story(input_dir, story_group_id);

%% Stack Stories
story_groups = [1 2 3 3 3 4];
fn_stack_stories(input_dir, story_ht, story_wt, story_groups, model_id)


