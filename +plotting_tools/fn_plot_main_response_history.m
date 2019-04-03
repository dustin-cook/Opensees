function [ ] = fn_plot_main_response_history( plot_dir, read_dir_opensees, node, analysis, eq_analysis_timespace, eq, eq_dt, direction, record_edp )
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

% Define EQ timespace
eq_timespace = eq_dt:eq_dt:(eq_dt*length(eq));

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
    node_ground_id = node.id(node.x == ground_x & node.z == ground_z & node.y == 0 & ~strcmp(node.fix,'[000000]'));
    node_second_east_id = node.id(node.x == east_x & node.z == east_z & node.story == 2);
    node_roof_east_id = node.id(node.x == east_x & node.z == east_z & node.story == 6);
    node_second_center_id = node.id(node.x == center_x & node.z == center_z & node.story == 2);
    node_roof_center_id = node.id(node.x == center_x & node.z == center_z & node.story == 6);
    ground_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_ground_id) '.mat']);
    second_east_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_second_east_id) '.mat']);
    roof_east_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_east_id) '.mat']);
    second_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_second_center_id) '.mat']);
    roof_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_center_id) '.mat']);

    fn_plot_response_history( ground_id_TH.nd_TH.(disp_tag), eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Ground Displacement ' direction ' (in)'], 30, record_edp.disp_TH_ground.(direction)  )
    fn_plot_response_history( ground_id_TH.nd_TH.(accel_tag), eq_timespace, eq, eq_dt, rh_plot_dir, ['Ground Acceleration ' direction ' (g)'], 30, record_edp.accel_TH_ground.(direction) )
    fn_plot_response_history( second_center_id_TH.nd_TH.(disp_tag), eq_analysis_timespace, eq, eq_dt, rh_plot_dir, ['Second Floor Displacement Center ' direction ' (in)'], 30, record_edp.disp_TH_second.(direction) )
    fn_plot_response_history( second_center_id_TH.nd_TH.(accel_tag), eq_timespace, eq, eq_dt, rh_plot_dir, ['Second Floor Acceleration Center ' direction ' (g)'], 30, record_edp.accel_TH_second.(direction) )
    fn_plot_response_history( roof_center_id_TH.nd_TH.(disp_tag), eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement Center ' direction ' (in)'], 30, record_edp.disp_TH_roof.(direction) )
    fn_plot_response_history( roof_center_id_TH.nd_TH.(disp_tag)/roof_ht, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Drift Center ' direction ' (in)'], 30, record_edp.disp_TH_roof.(direction)/roof_ht )
    fn_plot_response_history( roof_center_id_TH.nd_TH.(accel_tag), eq_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Acceleration Center ' direction ' (g)'], 30, record_edp.accel_TH_roof.(direction))
    if strcmp(direction,'z')
        fn_plot_response_history( roof_east_id_TH.nd_TH.(disp_tag), eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement East ' direction ' (in)'], 30, record_edp.disp_TH_roof_east.(direction) )
        fn_plot_response_history( roof_east_id_TH.nd_TH.(accel_tag), eq_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Acceleration East ' direction ' (g)'], 30, record_edp.accel_TH_roof_east.(direction) )
        % Torsional Response Roof
        relative_torsional_displacement = roof_east_id_TH.nd_TH.disp_z_TH - roof_center_id_TH.nd_TH.disp_z_TH;
        record_torsional_displacement = record_edp.disp_TH_roof_east.z - record_edp.disp_TH_roof.z;
        fn_plot_response_history( relative_torsional_displacement, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Twist (in)'], 30, record_torsional_displacement)
        % Torsional Response 1st story
        relative_torsional_displacement = second_east_id_TH.nd_TH.disp_z_TH - second_center_id_TH.nd_TH.disp_z_TH;
        record_torsional_displacement = record_edp.disp_TH_second_east.z - record_edp.disp_TH_second.z;
        fn_plot_response_history( relative_torsional_displacement, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['1st Story Twist (in)'], 30, record_torsional_displacement)
        % Torsional Response 1st story
        relative_torsional_displacement = (second_east_id_TH.nd_TH.disp_z_TH - second_center_id_TH.nd_TH.disp_z_TH)/174;
        record_torsional_displacement = (record_edp.disp_TH_second_east.z - record_edp.disp_TH_second.z)/174;
        fn_plot_response_history( relative_torsional_displacement, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['1st Story Twist Drift (rad)'], 30, record_torsional_displacement)
    end
else
    roof_nodes = node(node.y == roof_ht,:);
    mid_x = (max(roof_nodes.x) - min(roof_nodes.x)) / 2;
    mid_z = (max(roof_nodes.z) - min(roof_nodes.z)) / 2;
    for i = 1:height(roof_nodes)
        hyp_dist(i) = sqrt((roof_nodes.x(i) - mid_x)^2 + (roof_nodes.z(i) - mid_z)^2);
    end
    [~, idx] = min(hyp_dist);
    node_roof_center_id = roof_nodes.id(idx);
    roof_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_center_id) '.mat']);
    fn_plot_response_history( roof_center_id_TH.nd_TH.(disp_tag), eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Displacement Center ' direction ' (in)'], 30 )
    fn_plot_response_history( roof_center_id_TH.nd_TH.(disp_tag)/roof_ht, eq_analysis_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Drift Center ' direction], 30 )
    fn_plot_response_history( roof_center_id_TH.nd_TH.(accel_tag), eq_timespace, eq, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, ['Roof Acceleration Center ' direction ' (g)'], 30 )
end

%% Node time history comparison check (set up as sperate checker function
% node_roof = node(node.story == 3 & node.record_accel == 1,:);
% node_center = node_roof(node_roof.x == 671 & node_roof.z == 300,:);
% node_extreeme = node_roof(node_roof.x == 1642 & node_roof.z == 450,:);
% center_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_center.id) '.mat']);
% extreeme_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_extreeme.id) '.mat']);
% 
% hold on
% plot(eq_analysis_timespace,movmean(center_TH.nd_TH.accel_z_abs_TH,50))
% plot(eq_analysis_timespace,movmean(extreeme_TH.nd_TH.accel_z_abs_TH,50))
% xlim([0,15])
end

