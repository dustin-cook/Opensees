function [ ] = main_ASCE_41( analysis )
% Description: Main script facilitating an ASCE 41 teir 3 assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
import build_model.main_build_model

%% Begin Method
% Pull in database of available models
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);

% Select Model for the analysis 
model = model_table(model_table.id == analysis.model_id,:);

% Create Outputs Directory
output_dir = ['outputs/' model.name{1} filesep 'model data'];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% Build Model
main_build_model( model, analysis, output_dir )

end

