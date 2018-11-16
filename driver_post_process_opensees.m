% Post Process Opensees Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 12;
analysis.gm_id = 8;
analysis.name = 'linear';
analysis.nonlinear = 0;
analysis.type = 1;
analysis.full_recorders = 0;
analysis.run_eigen = 0;

%% Import Packages
import opensees.main_post_process_opensees

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable(['outputs' filesep model.name{1} filesep 'model data' filesep 'element.csv'],'ReadVariableNames',true);
node = readtable(['outputs' filesep model.name{1} filesep 'model data' filesep 'node.csv'],'ReadVariableNames',true);
story = readtable(['outputs' filesep model.name{1} filesep 'model data' filesep 'story.csv'],'ReadVariableNames',true);
hinge = [];
if analysis.nonlinear ~= 0
    hinge = readtable(['outputs' filesep model.name{1} filesep 'model data' filesep 'hinge.csv'],'ReadVariableNames',true);
end

dirs = {'x','y','z'};
for i = 1:length(dirs)
    if ground_motion_seq.(['eq_id_',dirs{i}]) ~= 0
        ground_motion.(dirs{i}) = gm_table(gm_table.id == ground_motion_seq.(['eq_id_',dirs{i}]),:);
    end
end

%% Run Post Processor
main_post_process_opensees( analysis, model, story, node, element, hinge, ground_motion, output_dir )
