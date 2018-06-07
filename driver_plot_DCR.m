% Combine DCR from load combos and plot
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = {'09DL', '11DL11LL'};

%% Import Packages
import asce_41.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
plot_dir = ['outputs' filesep model.name{1} filesep 'plots'];

for i = 1:length(analysis.name)
    output_dir = ['outputs' filesep model.name{1} filesep analysis.name{i}];
    element{i} = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
    node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);
end

%% Find DCR Envelope
element_combo = element{1};
for i = 1:length(analysis.name)
    element_combo.DCR_M = max([element_combo.DCR_M,element{i}.DCR_M],[],2);
    element_combo.DCR_V = max([element_combo.DCR_V,element{i}.DCR_V],[],2);
    element_combo.DCR_P = max([element_combo.DCR_P,element{i}.DCR_P],[],2);
end

%% Plot DCR
% moment
fn_plot_building( element_combo.DCR_M, element_combo, node, 'DCR_view_moment', plot_dir )

% shear
fn_plot_building( element_combo.DCR_V, element_combo, node, 'DCR_view_shear', plot_dir )

% axial
fn_plot_building( element_combo.DCR_P, element_combo, node, 'DCR_view_axial', plot_dir )