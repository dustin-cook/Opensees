clear all
close all
clc

% Basic Inputs
analysis.proceedure = 'NDP';
analysis.model_id = 11;

% Import Packages
import asce_41.main_hinge_properties

% Load Model Table
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Define Read and Write Directories
analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure];
read_dir = [analysis.out_dir filesep 'asce_41_data'];

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'joint_analysis.mat'])
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% Run Hinge Properties
[ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );