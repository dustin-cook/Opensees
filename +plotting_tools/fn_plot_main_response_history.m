function [ ] = fn_plot_main_response_history( plot_dir, read_dir_opensees, node, analysis, eq_analysis_timespace, eq, eq_dt, record_edp )
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
eq_timespace = eq_dt:eq_dt:(eq_dt*length(eq.x));
max_time2plot = 30;

%% Begin Method
roof_ht = max(node.y);
second_floor_ht = max(node.y(node.story == 1 & node.on_slab == 1));

% Plot recorded vs analysis signals
if ~isempty(record_edp)
    % Define Recored Drift
    if isfield(eq,'z')
        record_edp.drift_TH_roof.z = record_edp.disp_TH_roof.z/roof_ht;
        record_edp.drift_TH_second.z = record_edp.disp_TH_second.z/second_floor_ht;
    end
    record_edp.drift_TH_roof.x = record_edp.disp_TH_roof.x/roof_ht;
    record_edp.drift_TH_second.x = record_edp.disp_TH_second.x/second_floor_ht;
    

    % Define Accelerometer Nodes (Currently Specific to ICSB, NEED TO UPDATE)
    center_x = 671;
    center_z = 300;
    ground_x = 671;
    ground_z = 450;
    east_x = 1571;
    east_z = 300;
    node_ground_id = node.id(node.x == ground_x & node.z == ground_z & node.y == 0 & node.on_slab == 1);
    node_second_east_id = node.id(node.x == east_x & node.z == east_z & node.story == 2);
    node_roof_east_id = node.id(node.x == east_x & node.z == east_z & node.story == 6);
    node_second_center_id = node.id(node.x == center_x & node.z == center_z & node.story == 1);
    node_roof_center_id = node.id(node.x == center_x & node.z == center_z & node.story == 6);
    ground_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_ground_id) '.mat']);
    second_east_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_second_east_id) '.mat']);
    roof_east_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_east_id) '.mat']);
    second_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_second_center_id) '.mat']);
    roof_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_center_id) '.mat']);

    % Plot response histories
    if isfield(eq,'z')
        fn_plot_response_history( ground_id_TH.nd_TH.disp_x_TH, ground_id_TH.nd_TH.disp_z_TH, eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Ground Displacement (in)', max_time2plot, record_edp.disp_TH_ground )
        fn_plot_response_history( ground_id_TH.nd_TH.accel_x_abs_TH, ground_id_TH.nd_TH.accel_z_abs_TH, eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Ground Acceleration (g)', max_time2plot, record_edp.accel_TH_ground )
        fn_plot_response_history( second_center_id_TH.nd_TH.disp_x_TH, second_center_id_TH.nd_TH.disp_z_TH, eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Displacement Center (in)', max_time2plot, record_edp.disp_TH_second )
        fn_plot_response_history( second_center_id_TH.nd_TH.disp_x_TH/second_floor_ht, second_center_id_TH.nd_TH.disp_z_TH/second_floor_ht, eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Drift Center (in)', max_time2plot, record_edp.drift_TH_second )
        fn_plot_response_history( second_center_id_TH.nd_TH.accel_x_abs_TH, second_center_id_TH.nd_TH.accel_z_abs_TH, eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Acceleration Center (g)', max_time2plot, record_edp.accel_TH_second )
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH, roof_center_id_TH.nd_TH.disp_z_TH, eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Displacement Center (in)', max_time2plot, record_edp.disp_TH_roof )
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH/roof_ht, roof_center_id_TH.nd_TH.disp_z_TH/roof_ht, eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Drift Center (in)', max_time2plot, record_edp.drift_TH_roof )
        fn_plot_response_history( roof_center_id_TH.nd_TH.accel_x_abs_TH, roof_center_id_TH.nd_TH.accel_z_abs_TH, eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Acceleration Center (g)', max_time2plot, record_edp.accel_TH_roof)

        % East side z roof response
        fn_plot_response_history( roof_east_id_TH.nd_TH.disp_z_TH, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Displacement East (in)', max_time2plot, record_edp.disp_TH_roof_east.z )
        fn_plot_response_history( roof_east_id_TH.nd_TH.accel_z_abs_TH, [], eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Acceleration East (g)', max_time2plot, record_edp.accel_TH_roof_east.z )
        % East side z first story response
        fn_plot_response_history( second_east_id_TH.nd_TH.disp_z_TH, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Displacement East (in)', max_time2plot, record_edp.disp_TH_second_east.z )
        fn_plot_response_history( second_east_id_TH.nd_TH.accel_z_abs_TH, [], eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Acceleration East (g)', max_time2plot, record_edp.accel_TH_second_east.z )
        % Torsional Response Roof
        relative_torsional_displacement = roof_east_id_TH.nd_TH.disp_z_TH - roof_center_id_TH.nd_TH.disp_z_TH;
        record_torsional_displacement = record_edp.disp_TH_roof_east.z - record_edp.disp_TH_roof.z;
        fn_plot_response_history( relative_torsional_displacement, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Twist (in)', max_time2plot, record_torsional_displacement)
        % Torsional Response 1st story
        relative_torsional_displacement = second_east_id_TH.nd_TH.disp_z_TH - second_center_id_TH.nd_TH.disp_z_TH;
        record_torsional_displacement = record_edp.disp_TH_second_east.z - record_edp.disp_TH_second.z;
        fn_plot_response_history( relative_torsional_displacement, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, ['1st Story Twist (in)'], max_time2plot, record_torsional_displacement)
        % Torsional Response 1st story
        relative_torsional_displacement = (second_east_id_TH.nd_TH.disp_z_TH - second_center_id_TH.nd_TH.disp_z_TH)/second_floor_ht;
        record_torsional_displacement = (record_edp.disp_TH_second_east.z - record_edp.disp_TH_second.z)/second_floor_ht;
        fn_plot_response_history( relative_torsional_displacement, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, ['1st Story Twist Drift (rad)'], max_time2plot, record_torsional_displacement)
    else
        fn_plot_response_history( ground_id_TH.nd_TH.disp_x_TH, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Ground Displacement (in)', max_time2plot, record_edp.disp_TH_ground.x )
        fn_plot_response_history( ground_id_TH.nd_TH.accel_x_abs_TH, [], eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Ground Acceleration (g)', max_time2plot, record_edp.accel_TH_ground.x )
        fn_plot_response_history( second_center_id_TH.nd_TH.disp_x_TH, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Displacement Center (in)', max_time2plot, record_edp.disp_TH_second.x )
        fn_plot_response_history( second_center_id_TH.nd_TH.disp_x_TH/second_floor_ht, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Drift Center (in)', max_time2plot, record_edp.drift_TH_second.x )
        fn_plot_response_history( second_center_id_TH.nd_TH.accel_x_abs_TH, [], eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Second Floor Acceleration Center (g)', max_time2plot, record_edp.accel_TH_second.x )
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Displacement Center (in)', max_time2plot, record_edp.disp_TH_roof.x )
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH/roof_ht, [], eq_analysis_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Drift Center (in)', max_time2plot, record_edp.drift_TH_roof.x )
        fn_plot_response_history( roof_center_id_TH.nd_TH.accel_x_abs_TH, [], eq_timespace, eq.x, eq_dt, rh_plot_dir, 'Roof Acceleration Center (g)', max_time2plot, record_edp.accel_TH_roof.x )
    end
else
    node_roof_center_id = node.id(node.y == roof_ht & node.primary_story == 1);
    roof_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_center_id) '.mat']);
    if isfield(roof_center_id_TH.nd_TH,'disp_z_TH')
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH, roof_center_id_TH.nd_TH.disp_z_TH, eq_analysis_timespace, eq.x, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, 'Roof Displacement Center (in)', max_time2plot)
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH/roof_ht, roof_center_id_TH.nd_TH.disp_z_TH/roof_ht, eq_analysis_timespace, eq.x, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, 'Roof Drift Center (in)', max_time2plot)
    else
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH, [], eq_analysis_timespace, eq.x, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, 'Roof Displacement Center (in)', max_time2plot)
        fn_plot_response_history( roof_center_id_TH.nd_TH.disp_x_TH/roof_ht, [], eq_analysis_timespace, eq.x, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, 'Roof Drift Center (in)', max_time2plot)
    end
    if ~isempty(roof_center_id_TH.nd_TH.accel_x_abs_TH)
        fn_plot_response_history( roof_center_id_TH.nd_TH.accel_x_abs_TH, roof_center_id_TH.nd_TH.accel_z_abs_TH, eq_timespace, eq.x, eq_dt/analysis.initial_timestep_factor^2, rh_plot_dir, 'Roof Acceleration Center (g)', max_time2plot)
    end
end

end

