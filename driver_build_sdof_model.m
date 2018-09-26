% Build SDOF Model Tables
clear
close
clc

%% DEFINE INPTUTS
analysis.model_id = 10;

%% Initial Setup
import build_model.fn_build_sdof

% Load basic model info
model_table = readtable(['inputs' filesep 'sdof_models.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
% Create Outputs Directory
output_dir = ['outputs/' model.name{1} filesep 'model data'];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Build Model
fn_build_sdof( model, output_dir )