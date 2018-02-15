%% INPUTS
model_name = 'ICBS_model_3D';
story_ht = [174];
story_wt = [2200];
bay_length.x = [300 300 300 300 300];
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
col_id = [6 0];
beam_id = [8 0];
wall_id = [0 3];
grid_lines{1}.bays = 1:5;
grid_lines{2}.bays = 1;
direction = {'x' 'z'};

fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines)

%% Asseble Stories
story_group_id{1}.grid_id = [1 1 1 1 2 2 2 2];
story_group_id{1}.start_prim = [0 0 0 0 300 300 300 300];
story_group_id{1}.start_alt = [0 300 600 900 0 600 900 1200];
story_group_id{1}.direction = [1 1 1 1 3 3 3 3];

fn_assemble_story(input_dir, story_group_id);

%% Stack Stories
story_groups = [1];
fn_stack_stories(input_dir, story_ht, story_wt, story_groups, model_id)


