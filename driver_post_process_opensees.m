% Post Process Opensees Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 6;
analysis.name = 'test';
analysis.nonlinear = 0;

%% Import Packages
import opensees.main_post_process_opensees

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);
story = readtable([output_dir filesep 'story.csv'],'ReadVariableNames',true);

dirs = {'x','y','z'};
for i = 1:length(dirs)
    if ground_motion_seq.(['eq_id_',dirs{i}]) ~= 0
        ground_motion.(dirs{i}) = gm_table(gm_table.id == ground_motion_seq.(['eq_id_',dirs{i}]),:);
    end
end

%% Run Post Processor
main_post_process_opensees( analysis, model, story, node, element, ground_motion, output_dir )
