% Plot Analysis Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = 'NL_10DL10LL';
analysis.nonlinear = 1;

%% Import Packages
import plotting_tools.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
plot_dir = [output_dir filesep 'plots'];
load([output_dir filesep 'ASCE_data']);
% Max channel recordings
load([pwd filesep 'ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat'])

% Linear Procedures
if analysis.nonlinear == 0
    %% Plot DCR
    % envelope
    fn_plot_building( element.DCR_total, element, node, 'DCR_view_envelope', plot_dir )
    % moment
    fn_plot_building( element.DCR_M, element, node, 'DCR_view_moment', plot_dir )
    % shear
    fn_plot_building( element.DCR_V, element, node, 'DCR_view_shear', plot_dir )
    % axial
    fn_plot_building( element.DCR_P, element, node, 'DCR_view_axial', plot_dir )
end

for i = 1:length(dirs_ran)
    %% Plot EDP Profiles
    fn_plot_profile( [max(abs(eq.(dirs_ran(i)))); story.(['max_accel_' dirs_ran(i)])], [0;story.id], plot_dir, ['Acceleration Profile ' dirs_ran(i)], 'PFA (g)', 0.5, record_edp.max_accel.(dirs_ran(i)))
    fn_plot_profile( [0; story.(['max_disp_' dirs_ran(i)])], [0;story.id], plot_dir, ['Displacement Profile ' dirs_ran(i)], 'Displacement (in)', 10, record_edp.max_disp.(dirs_ran(i)) )
    fn_plot_profile( story.(['max_drift_' dirs_ran(i)]), story.id, plot_dir, ['Drift Profile ' dirs_ran(i)], 'IDR', 0.01 )

    % Plot Roof Response History
    fn_plot_response_history( node.(['disp_' dirs_ran(i) '_TH']), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Displacemnet ' dirs_ran(i) ' (in)'] )
    fn_plot_response_history( node.(['accel_' dirs_ran(i) '_abs_TH']), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Acceleration ' dirs_ran(i) ' (g)'] )
end
