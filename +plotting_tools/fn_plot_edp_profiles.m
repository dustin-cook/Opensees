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

%% Begin Method
if exist('record_edp','var')
    % Acceleration
    fn_plot_profile( [ground_motion.x.pga; story.max_accel_x], [ground_motion.z.pga; story.max_accel_z], [min(story.id)-1;story.id], edp_plot_dir, 'Acceleration Profile', 'Acceleration (g)', NaN, record_edp.max_accel)
    fn_plot_profile( [ground_motion.x.pga; story.max_accel_center_x], [ground_motion.z.pga; story.max_accel_center_z], [min(story.id)-1;story.id], edp_plot_dir, 'Acceleration Profile Center', 'Acceleration (g)', NaN, record_edp.max_accel_center)

    % Displacement
    fn_plot_profile( [0; story.max_disp_x], [0; story.max_disp_z], [min(story.id)-1;story.id], edp_plot_dir, 'Displacement Profile', 'Displacement (in)', target_disp_in, record_edp.max_disp)
    fn_plot_profile( [0; story.max_disp_center_x], [0; story.max_disp_center_z], [min(story.id)-1;story.id], edp_plot_dir, 'Displacement Profile Center', 'Displacement (in)', target_disp_in, record_edp.max_disp_center)
    
%     % Relative Displacement
%     recorded_roof_disp = record_edp.max_disp.(direction)(end);
%     analysis_roof_disp = story.(['max_disp_' direction])(end);
%     fn_plot_profile( [0; story.(['max_disp_' direction])]/analysis_roof_disp , [0;story.id], edp_plot_dir, ['Normalized Displacement Profile ' direction], 'Normalized Displacement', 1, record_edp.max_disp.(direction)/recorded_roof_disp  )
else 
    % Acceleration
    fn_plot_profile( [ground_motion.x.pga; story.max_accel_x], [ground_motion.z.pga; story.max_accel_z], [0;story.id], edp_plot_dir, 'Acceleration Profile', 'Peak Floor Acceleration (g)', NaN)
    fn_plot_profile( [ground_motion.x.pga; story.max_accel_center_x], [ground_motion.z.pga; story.max_accel_center_z], [0;story.id], edp_plot_dir, 'Acceleration Profile Center', 'Peak Floor Acceleration (g)', NaN)

    % Displacement
    fn_plot_profile( [0; story.max_disp_x], [0; story.max_disp_z], [0;story.id], edp_plot_dir, 'Displacement Profile', 'Peak Floor Displacement (in)', target_disp_in)
    fn_plot_profile( [0; story.max_disp_center_x], [0; story.max_disp_center_z], [0;story.id], edp_plot_dir, 'Displacement Profile Center', 'Peak Floor Displacement (in)', target_disp_in)
end
% Drift
fn_plot_profile( story.max_drift_x, story.max_drift_z, story.id, edp_plot_dir, 'Drift Profile', 'Interstory Drift', NaN )

end

