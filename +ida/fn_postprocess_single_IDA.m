function [] = fn_postprocess_single_IDA(analysis, model, gm_set_table)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Initial Setup
% Import Packages
import asce_41.*
import plotting_tools.*

% IDA directory
outputs_read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'scale_' num2str(analysis.single_gm_scale) '/' 'GM_' num2str(analysis.single_gm_set_id) '_' num2str(analysis.single_gm_pair_id)];
summary_read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(analysis.single_gm_scale) '/' 'GM_' num2str(analysis.single_gm_set_id) '_' num2str(analysis.single_gm_pair_id)];
plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'IDA Plots' '/' 'scale_' num2str(analysis.single_gm_scale) '/' 'GM_' num2str(analysis.single_gm_set_id) '_' num2str(analysis.single_gm_pair_id)];

%% Load in Data
% Pull in Element Database
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% Load analysis data
load([outputs_read_dir filesep 'element_analysis.mat'])
load([summary_read_dir filesep 'hinge_analysis.mat'])
load([outputs_read_dir filesep 'node_analysis.mat'])
load([summary_read_dir filesep 'story_analysis.mat'])
converge_tol_data = readmatrix([outputs_read_dir filesep 'converge_tol_file.txt']);

%% Define Variables and Post Process Hinges
% Define this ground motion
gm_set = gm_set_table(gm_set_table.set_id == analysis.single_gm_set_id,:);
if analysis.single_gm_pair_id == 1
    ground_motion.x.pga = gm_set.pga(gm_set.pair == 1);
    ground_motion.z.pga = gm_set.pga(gm_set.pair == 2);
elseif analysis.single_gm_pair_id == 2
    ground_motion.x.pga = gm_set.pga(gm_set.pair == 2);
    ground_motion.z.pga = gm_set.pga(gm_set.pair == 1);
end
ground_motion.x.eq_dt = gm_set.eq_dt(1);
eq_analysis_timespace = converge_tol_data(:,1);
eq.x = gm_set.eq_dt(1):gm_set.eq_dt(1):(gm_set.eq_dt(1)*gm_set.eq_length(1));

% filter OOP hinge out of table
hinge = hinge(~strcmp(hinge.direction,'oop'),:);

% calculate hinge demand to capacity ratios
[ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge, outputs_read_dir, node );

%% Plot Results for this IDA gm
% Run Plotters
fn_plot_element_scatter( element, 'column', story, hinge, plot_dir )
fn_plot_element_scatter( element, 'beam', story, hinge, plot_dir )
fn_plot_element_scatter( element, 'wall', story, hinge, plot_dir )
fn_plot_edp_profiles( plot_dir, ground_motion, model, story, NaN)
fn_plot_main_response_history( plot_dir, outputs_read_dir, node, analysis, eq_analysis_timespace, eq, ground_motion.x.eq_dt )
analysis.hinge_stories_2_plot = 1;
fn_plot_hinge_response( plot_dir, outputs_read_dir, hinge, element, ele_prop_table, node, analysis.hinge_stories_2_plot, eq_analysis_timespace )

end

