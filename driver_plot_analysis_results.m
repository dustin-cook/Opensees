% Plot Analysis Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '11DL11LL';
analysis.nonlinear = 0;

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
    fn_plot_building( element.DCR_total, element, node, 'DCR_view_envelope_ext', plot_dir, '3D', 'linear', 'ext' )
    % moment
    fn_plot_building( element.DCR_M, element, node, 'DCR_view_moment_ext', plot_dir, '3D', 'linear', 'ext' )
    % shear
    fn_plot_building( element.DCR_V, element, node, 'DCR_view_shear_ext', plot_dir, '3D', 'linear', 'ext' )
    % axial
    fn_plot_building( element.DCR_P, element, node, 'DCR_view_axial_ext', plot_dir, '3D', 'linear', 'ext' )
    
    %% Plot DCR raw
    % envelope
    fn_plot_building( element.DCR_total_raw, element, node, 'DCR_view_envelope_ext_raw', plot_dir, '3D', 'raw', 'ext' )
    % moment
    fn_plot_building( element.DCR_M_raw, element, node, 'DCR_view_moment_ext_raw', plot_dir, '3D', 'raw', 'ext' )
    % shear
    fn_plot_building( element.DCR_V_raw, element, node, 'DCR_view_shear_ext_raw', plot_dir, '3D', 'raw', 'ext' )
    % axial
    fn_plot_building( element.DCR_P_raw, element, node, 'DCR_view_axial_ext_raw', plot_dir, '3D', 'raw', 'ext' )
    
    %% Plot DCR
    % envelope
    fn_plot_building( element.DCR_total, element, node, 'DCR_view_envelope_int', plot_dir, '3D', 'linear', 'int' )
    % moment
    fn_plot_building( element.DCR_M, element, node, 'DCR_view_moment_int', plot_dir, '3D', 'linear', 'int' )
    % shear
    fn_plot_building( element.DCR_V, element, node, 'DCR_view_shear_int', plot_dir, '3D', 'linear', 'int' )
    % axial
    fn_plot_building( element.DCR_P, element, node, 'DCR_view_axial_int', plot_dir, '3D', 'linear', 'int' )
    
    %% Plot DCR raw
    % envelope
    fn_plot_building( element.DCR_total_raw, element, node, 'DCR_view_envelope_int_raw', plot_dir, '3D', 'raw', 'int' )
    % moment
    fn_plot_building( element.DCR_M_raw, element, node, 'DCR_view_moment_int_raw', plot_dir, '3D', 'raw', 'int' )
    % shear
    fn_plot_building( element.DCR_V_raw, element, node, 'DCR_view_shear_int_raw', plot_dir, '3D', 'raw', 'int' )
    % axial
    fn_plot_building( element.DCR_P_raw, element, node, 'DCR_view_axial_int_raw', plot_dir, '3D', 'raw', 'int' )
else
    %% Plot accpetance
    fn_plot_building_nl( hinge, element, node, 'Acceptance Plot', plot_dir)
end

for i = 1:length(dirs_ran)
    %% Plot EDP Profiles
    fn_plot_profile( [max(abs(eq.(dirs_ran(i)))); story.(['max_accel_' dirs_ran(i)])], [0;story.id], plot_dir, ['Acceleration Profile ' dirs_ran(i)], 'PFA (g)', 0.5, record_edp.max_accel.(dirs_ran(i)))
    fn_plot_profile( [0; story.(['max_disp_' dirs_ran(i)])], [0;story.id], plot_dir, ['Displacement Profile ' dirs_ran(i)], 'Displacement (in)', 10, record_edp.max_disp.(dirs_ran(i)) )
    fn_plot_profile( story.(['max_drift_' dirs_ran(i)]), story.id, plot_dir, ['Drift Profile ' dirs_ran(i)], 'SDR', 0.01 )
    fn_plot_profile( [0; story.(['max_disp_' dirs_ran(i) '_ASCE'])], [0;story.id], plot_dir, ['ASCE Displacement Profile ' dirs_ran(i)], 'Displacement (in)', 10, record_edp.max_disp.(dirs_ran(i)) )
    fn_plot_profile( story.(['max_drift_' dirs_ran(i) '_ASCE']), story.id, plot_dir, ['ASCE Drift Profile ' dirs_ran(i)], 'SDR', 0.01 )
    
    % Plot specific TH comparisons
    if strcmp(model.dimension,'3D')
        node_ground = node(node.x == 1200 & node.y == 0 & node.z == 300,:);
        node_second_east = node(node.x == 1500 & node.y == 174 & node.z == 300,:);
        node_roof_east = node(node.x == 1500 & node.y == 822 & node.z == 300,:);
        node_second_center = node(node.x == 900 & node.y == 174 & node.z == 300,:);
        node_roof_center = node(node.x == 900 & node.y == 822 & node.z == 300,:);
        fn_plot_response_history( node_ground.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Displacemnet ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_ground.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Acceleration ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_second_east.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Displacemnet East ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_second_east.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Acceleration East ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_roof_east.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Displacemnet East ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_roof_east.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Acceleration East ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_second_center.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Displacemnet Center ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_second_center.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Acceleration Center ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_roof_center.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Displacemnet Center ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_roof_center.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Acceleration Center ' dirs_ran(i) ' (g)'] )
    elseif strcmp(model.dimension,'2D')
        node_ground = node(node.x == 1200 & node.y == 0 & node.id < 300,:);
        node_second_center = node(node.x == 900 & node.y == 174 & node.id < 300,:);
        node_roof_center = node(node.x == 900 & node.y == 822 & node.id < 300,:);
        fn_plot_response_history( node_ground.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Displacemnet ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_ground.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Acceleration ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_second_center.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Displacemnet Center ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_second_center.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Acceleration Center ' dirs_ran(i) ' (g)'] )
        fn_plot_response_history( node_roof_center.(['disp_' dirs_ran(i) '_TH'])-node_roof_center.(['disp_' dirs_ran(i) '_TH'])(1), record_edp.disp_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Displacemnet Center ' dirs_ran(i) ' (in)'] )
        fn_plot_response_history( node_roof_center.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Acceleration Center ' dirs_ran(i) ' (g)'] )
    end
end

