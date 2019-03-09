function [ ] = fn_plot_main_response_history( plot_dir, node, analysis, eq_analysis_timespace, eq, eq_dt, direction, record_edp )
% Description: Fn to plot response histories of the analysis

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import plotting_tools.fn_plot_response_history

% Define plot directory
rh_plot_dir = [plot_dir filesep 'Response History Plots'];

%% Begin Method
% Plot recorded vs analysis signals
disp_tag = ['disp_' direction '_TH'];
accel_tag = ['accel_' direction '_abs_TH'];
roof_ht = max(node.y);

if exist('record_edp','var')
    % Define Accelerometer Nodes (Currently Specific to ICSB, NEED TO UPDATE)
    center_x = 671;
    center_z = 300;
    ground_x = 1271;
    ground_z = 450;
    east_x = 1571;
    east_z = 300;
    node_ground = node(node.x == ground_x & node.z == ground_z & node.y == 0 & ~strcmp(node.fix,'[000000]'),:);
    node_second_east = node(node.x == east_x & node.z == east_z & node.story == 2,:);
    node_roof_east = node(node.x == east_x & node.z == east_z & node.story == 6,:);
    node_second_center = node(node.x == center_x & node.z == center_z & node.story == 2,:);
    node_roof_center = node(node.x == center_x & node.z == center_z & node.story == 6,:);

    fn_plot_response_history( node_ground.(disp_tag){1}, eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Ground Displacement ' direction ' (in)'], 30, record_edp.disp_TH_ground.(direction)  )
    fn_plot_response_history( node_ground.(accel_tag){1}, eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Ground Acceleration ' direction ' (g)'], 30, record_edp.accel_TH_ground.(direction) )
    fn_plot_response_history( node_second_center.(disp_tag){1}, eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Second Floor Displacement Center ' direction ' (in)'], 30, record_edp.disp_TH_second.(direction) )
    fn_plot_response_history( node_second_center.(accel_tag){1}, eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Second Floor Acceleration Center ' direction ' (g)'], 30, record_edp.accel_TH_second.(direction) )
    fn_plot_response_history( node_roof_center.(disp_tag){1}, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement Center ' direction ' (in)'], 30, record_edp.disp_TH_roof.(direction) )
    fn_plot_response_history( node_roof_center.(disp_tag){1}/roof_ht, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Drift Center ' direction ' (in)'], 30, record_edp.disp_TH_roof.(direction)/roof_ht )
    fn_plot_response_history( node_roof_center.(accel_tag){1}, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Acceleration Center ' direction ' (g)'], 30, record_edp.accel_TH_roof.(direction))
    if strcmp(direction,'z')
        fn_plot_response_history( node_roof_east.(disp_tag){1}, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement East ' direction ' (in)'], 30, record_edp.disp_TH_roof_east.(direction) )
        % Torsional Response
        relative_torsional_displacement = node_roof_east.disp_z_TH{1} - node_roof_center.disp_z_TH{1};
        record_torsional_displacement = record_edp.disp_TH_roof_east.z - record_edp.disp_TH_roof.z;
        fn_plot_response_history( relative_torsional_displacement, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Relative Torsional Displacement (in)'], 30, record_torsional_displacement)
    end
else
    roof_nodes = node(node.y == roof_ht,:);
    mid_x = (max(roof_nodes.x) - min(roof_nodes.x)) / 2;
    mid_z = (max(roof_nodes.z) - min(roof_nodes.z)) / 2;
    for i = 1:height(roof_nodes)
        hyp_dist(i) = sqrt((roof_nodes.x(i) - mid_x)^2 + (roof_nodes.z(i) - mid_z)^2);
    end
    [~, idx] = min(hyp_dist);
    node_roof_center = roof_nodes(idx,:);
    fn_plot_response_history( node_roof_center.(disp_tag){1}, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement Center ' direction ' (in)'], 30 )
    fn_plot_response_history( node_roof_center.(disp_tag){1}/roof_ht, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Drift Center ' direction], 30 )
    fn_plot_response_history( node_roof_center.(accel_tag){1}, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Acceleration Center ' direction ' (g)'], 30 )
end
end

