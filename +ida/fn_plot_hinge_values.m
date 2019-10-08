function [ ] = fn_plot_hinge_values(analysis, model, IDA_scale_factors, gm_set_table)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
import plotting_tools.*

plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'IDA Plots'];
if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end

%% Load data
% Hinges from IDA
outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(IDA_scale_factors(1)) '/' 'GM_' num2str(gm_set_table.set_id(1)) '_' num2str(gm_set_table.pair(1))];
hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
load(hinge_file)

% Nodes from analysis
node_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'asce_41_data'];
node_file = [node_dir filesep 'node_analysis.mat'];
load(node_file)

% Elements from analsysis
ele_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'asce_41_data'];
ele_file = [ele_dir filesep 'element_analysis.mat'];
load(ele_file)

% ICSB Demand Data
hin_icsb_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'asce_41_data'];
hin_icsb_file = [hin_icsb_dir filesep 'hinge_analysis.mat'];
hin_icsb = load(hin_icsb_file);

ele_sides = {'Base', 'Top'};
for i = 1:2
    %% Collect Hinge Data for Column Bases
    ele_2_use = element(element.story == 1 & strcmp(element.type,'column'),:);
    node_2_use = node(ismember(node.id,ele_2_use.node_1),:);
    hinge_2_use = hinge(ismember(hinge.element_id,ele_2_use.id) & hinge.ele_side == i & strcmp(hinge.direction,'primary'),:);
    hinge_icsb_2_use = hin_icsb.hinge(ismember(hin_icsb.hinge.element_id,ele_2_use.id) & hin_icsb.hinge.ele_side == i & strcmp(hin_icsb.hinge.direction,'primary'),:);
    hin_total_demands = hinge_icsb_2_use.b_ratio .* hinge_2_use.b_value_tot;

    %% Plot Plan view scatters
    % b values
    plot_name = ['Plan View - Column b values - ' ele_sides{i}];
    plot_txt = ['1st Story - ' ele_sides{i} ' of Columns'];
    x_lab = 'Total Rotation at b of ASCE 41';
    fn_plot_plan_scatter( node_2_use, hinge_2_use.b_value_tot, plot_dir, plot_name, 0, plot_txt, x_lab, [], 0.06 )

    % CP values
    plot_name = ['Plan View - Column CP values - ' ele_sides{i}];
    plot_txt = ['1st Story - ' ele_sides{i} ' of Columns'];
    x_lab = 'Total Rotation at CP of ASCE 41';
    fn_plot_plan_scatter( node_2_use, hinge_2_use.cp_value_tot, plot_dir, plot_name, 0, plot_txt, x_lab, [], 0.06 )

    % Eurocode NC values
    plot_name = ['Plan View - Column Euro NC values - ' ele_sides{i}];
    plot_txt = ['1st Story - ' ele_sides{i} ' of Columns'];
    x_lab = 'Total Rotation at NC of EuroCode';
    fn_plot_plan_scatter( node_2_use, hinge_2_use.euro_th_NC_value, plot_dir, plot_name, 0, plot_txt, x_lab, [], 0.06 )
    
    % Eurocode NC values
    plot_name = ['Plan View - Column Euro NC ratios - ' ele_sides{i}];
    plot_txt = ['1st Story - ' ele_sides{i} ' of Columns'];
    x_lab = 'Max(\theta) /NC_{Eurocode}';
    euro_ratios = hin_total_demands ./ hinge_2_use.euro_th_NC_value;
    fn_plot_plan_scatter( node_2_use, euro_ratios, plot_dir, plot_name, 0, plot_txt, x_lab, [], 2 )
end
end
