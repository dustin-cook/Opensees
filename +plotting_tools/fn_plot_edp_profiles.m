function [ ] = fn_plot_edp_profiles( plot_dir, ground_motion, model, story, target_disp_in, record_edp )
% Description: Fn to plot edp profiles of the analysis

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import plotting_tools.fn_plot_profile

% Plot Directory
edp_plot_dir = [plot_dir filesep 'EDP Profiles'];

%% Caclulate Profiles
if strcmp('max_accel_x',story.Properties.VariableNames)
    edp_accel_x = [ground_motion.x.pga; story.max_accel_x];
    edp_accel_x_center = [ground_motion.x.pga; story.max_accel_center_x];
end
edp_disp_x = [0; story.max_disp_x];
edp_disp_x_center = [0; story.max_disp_center_x];
if isfield(ground_motion,'z')
    if strcmp('max_accel_z',story.Properties.VariableNames)
        edp_accel_z = [ground_motion.z.pga; story.max_accel_z];
        edp_accel_z_center = [ground_motion.z.pga; story.max_accel_center_z];
    end
    edp_disp_z = [0; story.max_disp_z];
    edp_disp_z_center = [0; story.max_disp_center_z];
else
    edp_accel_z = [];
    edp_accel_z_center = [];
    edp_disp_z = [];
    edp_disp_z_center = [];
end

%% Plot Profiles
if exist('record_edp','var')
    % Acceleration
    if strcmp('max_accel_x',story.Properties.VariableNames)
        fn_plot_profile( edp_accel_x, edp_accel_z, [min(story.id)-1;story.id], edp_plot_dir, 'Acceleration Profile', 'Acceleration (g)', NaN, record_edp.max_accel)
        fn_plot_profile( edp_accel_x_center, edp_accel_z_center, [min(story.id)-1;story.id], edp_plot_dir, 'Acceleration Profile Center', 'Acceleration (g)', NaN, record_edp.max_accel_center)
    end
    % Displacement
    fn_plot_profile( edp_disp_x, edp_disp_z, [min(story.id)-1;story.id], edp_plot_dir, 'Displacement Profile', 'Displacement (in)', target_disp_in, record_edp.max_disp)
    if ~any(isnan(edp_disp_x_center))
        fn_plot_profile( edp_disp_x_center, edp_disp_z_center, [min(story.id)-1;story.id], edp_plot_dir, 'Displacement Profile Center', 'Displacement (in)', target_disp_in, record_edp.max_disp_center)
    end
else 
    % Acceleration
    if strcmp('max_accel_x',story.Properties.VariableNames)
        fn_plot_profile( edp_accel_x, edp_accel_z, [0;story.id], edp_plot_dir, 'Acceleration Profile', 'Peak Floor Acceleration (g)', NaN)
        fn_plot_profile( edp_accel_x_center, edp_accel_z_center, [0;story.id], edp_plot_dir, 'Acceleration Profile Center', 'Peak Floor Acceleration (g)', NaN)
    end
    % Displacement
    fn_plot_profile( edp_disp_x, edp_disp_z, [0;story.id], edp_plot_dir, 'Displacement Profile', 'Peak Floor Displacement (in)', target_disp_in)
    if ~any(isnan(edp_disp_x_center))
        fn_plot_profile( edp_disp_x_center, edp_disp_z_center, [0;story.id], edp_plot_dir, 'Displacement Profile Center', 'Peak Floor Displacement (in)', target_disp_in)
    end
end
% Drift
if isfield(ground_motion,'z')
    fn_plot_profile( story.max_drift_x, story.max_drift_z, story.id, edp_plot_dir, 'Drift Profile', 'Interstory Drift', NaN )
else
    fn_plot_profile( story.max_drift_x, [], story.id, edp_plot_dir, 'Drift Profile', 'Interstory Drift', NaN )
end

end

