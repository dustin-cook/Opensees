function [ ] = fn_plot_edp_profiles( plot_dir, pga, model, story, target_disp_in, direction, record_edp )
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
    fn_plot_profile( [pga; story.(['max_accel_' direction])], [0;story.id], edp_plot_dir, ['Acceleration Profile ' direction], 'PFA (g)', 0.8, NaN, model.num_stories, record_edp.max_accel.(direction))
    fn_plot_profile( [pga; story.(['max_accel_center_' direction])], [0;story.id], edp_plot_dir, ['Acceleration Profile ' direction ' Center'], 'PFA (g)', 0.8, NaN, model.num_stories, record_edp.max_accel_center.(direction))

    % Displacement
    fn_plot_profile( [0; story.(['max_disp_' direction])], [0;story.id], edp_plot_dir, ['Displacement Profile ' direction], 'Displacement (in)', 10, target_disp_in, model.num_stories, record_edp.max_disp.(direction))
    fn_plot_profile( [0; story.(['max_disp_center_' direction])], [0;story.id], edp_plot_dir, ['Displacement Profile ' direction ' Center'], 'Displacement (in)', 10, target_disp_in, model.num_stories, record_edp.max_disp_center.(direction))
    
    % Relative Displacement
    recorded_roof_disp = record_edp.max_disp.(direction)(end);
    analysis_roof_disp = story.(['max_disp_' direction])(end);
    fn_plot_profile( [0; story.(['max_disp_' direction])]/analysis_roof_disp , [0;story.id], edp_plot_dir, ['Normalized Displacement Profile ' direction], 'Normalized Displacement', 1, record_edp.max_disp.(direction)/recorded_roof_disp  )
else 
    % Acceleration
    fn_plot_profile( [pga; story.(['max_accel_' direction])], [0;story.id], edp_plot_dir, ['Acceleration Profile ' direction], 'PFA (g)', 0.8, NaN, model.num_stories)
    fn_plot_profile( [pga; story.(['max_accel_center_' direction])], [0;story.id], edp_plot_dir, ['Acceleration Profile ' direction ' Center'], 'PFA (g)', 0.8, NaN, model.num_stories)

    % Displacement
    fn_plot_profile( [0; story.(['max_disp_' direction])], [0;story.id], edp_plot_dir, ['Displacement Profile ' direction], 'Displacement (in)', 10, target_disp_in, model.num_stories)
    fn_plot_profile( [0; story.(['max_disp_center_' direction])], [0;story.id], edp_plot_dir, ['Displacement Profile ' direction ' Center'], 'Displacement (in)', 10, target_disp_in, model.num_stories)
end
% Drift
fn_plot_profile( story.(['max_drift_' direction]), story.id, edp_plot_dir, ['Drift Profile ' direction], 'IDR', 0.05, NaN, model.num_stories )

end

