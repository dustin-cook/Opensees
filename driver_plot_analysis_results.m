% Plot Analysis Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 12;
analysis.gm_id = 8;
analysis.name = 'nonlinear';
analysis.nonlinear = 1;
analysis.type = 1;
analysis.pushover_direction = 'x';
analysis.initial_timestep_factor = 1;
plot_asce = 0;

%% Import Packages
import plotting_tools.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
plot_dir = [output_dir filesep 'plots'];
load([output_dir filesep 'node_analysis.mat'])
load([output_dir filesep 'story_analysis.mat'])
load([output_dir filesep 'gm_data.mat'])


%% Pushover Analysis
if analysis.type == 2
    control_nodes = node(node.primary_story == 1,:);
    base_nodes = node(node.y == 0,:);
    
    % Calulate base shear
    base_shear = abs(sum(base_nodes.(['reaction_' analysis.pushover_direction '_TH'])));
    
    % Calculate roof disp
    roof_node = control_nodes(control_nodes.y == max(control_nodes.y),:);
    roof_disp = roof_node.(['disp_' analysis.pushover_direction '_TH']);
    
    % Plot Roof Disp Pushover
    plot(roof_disp,base_shear/1000)
    ylabel('Total Base Shear (k)')
    xlabel('Roof Displacement (in)')
    plot_dir = [output_dir filesep 'Pushover_Plots'];
    plot_name = ['Roof Pushover - ' analysis.pushover_direction];
    fn_format_and_save_plot( plot_dir, plot_name, 2 )
    
    % Plot story Drift Pushover
    figure
    hold on
    for i = 1:height(control_nodes)
        story_node = control_nodes(i,:);
        story_disp(i,:) = story_node.(['disp_' analysis.pushover_direction '_TH']);
        if i == 1
            rel_story_disp = story_disp;
            story_drift = rel_story_disp ./ story_node.y;
        else
            rel_story_disp = story_disp(i,:) - story_disp(i-1,:);
            story_drift = rel_story_disp ./ (control_nodes.y(i) - control_nodes.y(i-1));
        end
        plot(story_drift,base_shear/1000,'DisplayName',['Story - ' num2str(i)'])
    end
    ylabel('Total Base Shear (k)')
    xlabel('Story Drift (in)')
    plot_dir = [output_dir filesep 'Pushover_Plots'];
    plot_name = ['Story Pushover - ' analysis.pushover_direction];
    fn_format_and_save_plot( plot_dir, plot_name, 1 )
    
%% Dynamic Analysis
elseif analysis.type == 1
    %% Load Data
    load([output_dir filesep 'element_TH.mat']);
    % Max channel recordings
    load([pwd filesep 'ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat'])
    
    if plot_asce
        %% Load Data
        load([output_dir filesep 'ASCE_data']);
        load([output_dir filesep 'element_PM.mat'])
    
        %% Linear Procedures
        if analysis.nonlinear == 0  
            if strcmp(model.dimension,'3D')% for 3D linear analysis
                %% Plot DCR
                % envelope
                fn_plot_building( element.DCR_max_all, element, node, 'DCR_view_envelope_ext', plot_dir, '3D', 'linear', 'ext' )
                % moment
                fn_plot_building( element.DCR_max_M, element, node, 'DCR_view_moment_ext', plot_dir, '3D', 'linear', 'ext' )
                % shear
                fn_plot_building( element.DCR_max_V, element, node, 'DCR_view_shear_ext', plot_dir, '3D', 'linear', 'ext' )
                % axial
                fn_plot_building( element.DCR_max_P, element, node, 'DCR_view_axial_ext', plot_dir, '3D', 'linear', 'ext' )

                %% Plot DCR raw
                % envelope
                fn_plot_building( element.DCR_raw_max_all, element, node, 'DCR_view_envelope_ext_raw', plot_dir, '3D', 'raw', 'ext' )
                % moment
                fn_plot_building( element.DCR_raw_max_M, element, node, 'DCR_view_moment_ext_raw', plot_dir, '3D', 'raw', 'ext' )
                % shear
                fn_plot_building( element.DCR_raw_max_V, element, node, 'DCR_view_shear_ext_raw', plot_dir, '3D', 'raw', 'ext' )
                % axial
                fn_plot_building( element.DCR_raw_max_P, element, node, 'DCR_view_axial_ext_raw', plot_dir, '3D', 'raw', 'ext' )

                %% Plot DCR
                % envelope
                fn_plot_building( element.DCR_max_all, element, node, 'DCR_view_envelope_int', plot_dir, '3D', 'linear', 'int' )
                % moment
                fn_plot_building( element.DCR_max_M, element, node, 'DCR_view_moment_int', plot_dir, '3D', 'linear', 'int' )
                % shear
                fn_plot_building( element.DCR_max_V, element, node, 'DCR_view_shear_int', plot_dir, '3D', 'linear', 'int' )
                % axial
                fn_plot_building( element.DCR_max_P, element, node, 'DCR_view_axial_int', plot_dir, '3D', 'linear', 'int' )

                %% Plot DCR raw
                % envelope
                fn_plot_building( element.DCR_raw_max_all, element, node, 'DCR_view_envelope_int_raw', plot_dir, '3D', 'raw', 'int' )
                % moment
                fn_plot_building( element.DCR_raw_max_M, element, node, 'DCR_view_moment_int_raw', plot_dir, '3D', 'raw', 'int' )
                % shear
                fn_plot_building( element.DCR_raw_max_V, element, node, 'DCR_view_shear_int_raw', plot_dir, '3D', 'raw', 'int' )
                % axial
                fn_plot_building( element.DCR_raw_max_P, element, node, 'DCR_view_axial_int_raw', plot_dir, '3D', 'raw', 'int' )
            else
                % envelope
                fn_plot_building_2D( element.DCR_max_all, element, node, 'DCR_view_envelope_ext', plot_dir,  'linear' )
                % moment
                fn_plot_building_2D( element.DCR_max_M, element, node, 'DCR_view_moment_ext', plot_dir, 'linear' )
                % shear
                fn_plot_building_2D( element.DCR_max_V, element, node, 'DCR_view_shear_ext', plot_dir, 'linear' )
                % axial
                fn_plot_building_2D( element.DCR_max_P, element, node, 'DCR_view_axial_ext', plot_dir, 'linear' )

                %% Plot DCR raw
                % envelope
                fn_plot_building_2D( element.DCR_raw_max_all, element, node, 'DCR_view_envelope_ext_raw', plot_dir,  'raw' )
                % moment
                fn_plot_building_2D( element.DCR_raw_max_M, element, node, 'DCR_view_moment_ext_raw', plot_dir, 'raw' )
                % shear
                fn_plot_building_2D( element.DCR_raw_max_V, element, node, 'DCR_view_shear_ext_raw', plot_dir,'raw')
                % axial
                fn_plot_building_2D( element.DCR_raw_max_P, element, node, 'DCR_view_axial_ext_raw', plot_dir, 'raw' )
            end
        else
            %% Plot accpetance
            fn_plot_building_nl( hinge, element, node, 'Acceptance Plot', plot_dir)
        end
    end
    
    for i = 1:length(dirs_ran)
        if ~strcmp(dirs_ran(i),'y') % Update way I am doing this directional thing
            %% Plot EDP Profiles
            fn_plot_profile( [max(abs(eq.(dirs_ran{i}))); story.(['max_accel_' dirs_ran{i}])], [0;story.id], plot_dir, ['Acceleration Profile ' dirs_ran{i}], 'PFA (g)', 0.5, record_edp.max_accel.(dirs_ran{i}))
            fn_plot_profile( [0; story.(['max_disp_' dirs_ran{i}])], [0;story.id], plot_dir, ['Displacement Profile ' dirs_ran{i}], 'Displacement (in)', 10, record_edp.max_disp.(dirs_ran{i}) )
            fn_plot_profile( story.(['max_drift_' dirs_ran{i}]), story.id, plot_dir, ['Drift Profile ' dirs_ran{i}], 'SDR', 0.01 )
            if plot_asce
                fn_plot_profile( [0; story.(['max_disp_' dirs_ran{i} '_ASCE'])], [0;story.id], plot_dir, ['ASCE Displacement Profile ' dirs_ran{i}], 'Displacement (in)', 10, record_edp.max_disp.(dirs_ran{i}) )
                fn_plot_profile( story.(['max_drift_' dirs_ran{i} '_ASCE']), story.id, plot_dir, ['ASCE Drift Profile ' dirs_ran{i}], 'SDR', 0.01 )
            end
            
            % Plot specific TH comparisons
            if strcmp(model.dimension,'3D')
%                 node_ground = node(node.primary_story == 1,:);
%                 node_second_east = node(node.primary_story == 2,:);
%                 node_roof_east = node(node.primary_story == 6,:);
%                 node_second_center = node(node.primary_story == 2,:);
                node_roof_center = node(node.primary_story == 1 & node.story == 3,:);
%                 fn_plot_response_history( node_ground.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Displacemnet ' dirs_ran(i) ' (in)'] )
%                 fn_plot_response_history( node_ground.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_ground.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Ground Acceleration ' dirs_ran(i) ' (g)'] )
%                 fn_plot_response_history( node_second_east.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Displacemnet East ' dirs_ran(i) ' (in)'] )
%                 fn_plot_response_history( node_second_east.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Acceleration East ' dirs_ran(i) ' (g)'] )
%                 fn_plot_response_history( node_roof_east.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Displacemnet East ' dirs_ran(i) ' (in)'] )
%                 fn_plot_response_history( node_roof_east.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_roof.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Roof Acceleration East ' dirs_ran(i) ' (g)'] )
%                 fn_plot_response_history( node_second_center.(['disp_' dirs_ran(i) '_TH']), record_edp.disp_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Displacemnet Center ' dirs_ran(i) ' (in)'] )
%                 fn_plot_response_history( node_second_center.(['accel_' dirs_ran(i) '_abs_TH']), record_edp.accel_TH_second.(dirs_ran(i)), eq.(dirs_ran(i)), ground_motion.(dirs_ran(i)).eq_dt, plot_dir, ['Second Floor Acceleration Center ' dirs_ran(i) ' (g)'] )
                fn_plot_response_history( node_roof_center.(['disp_' dirs_ran{i} '_TH']), record_edp.disp_TH_roof.(dirs_ran{i}), eq.(dirs_ran{i}), ground_motion.(dirs_ran{i}).eq_dt/analysis.initial_timestep_factor^2, plot_dir, ['Roof Displacemnet Center ' dirs_ran{i} ' (in)'] )
                fn_plot_response_history( node_roof_center.(['accel_' dirs_ran{i} '_abs_TH']), record_edp.accel_TH_roof.(dirs_ran{i}), eq.(dirs_ran{i}), ground_motion.(dirs_ran{i}).eq_dt/analysis.initial_timestep_factor^2, plot_dir, ['Roof Acceleration Center ' dirs_ran{i} ' (g)'] )
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
    end
    %% Plot Hinge Response
    if analysis.nonlinear ~= 0
        load([output_dir filesep 'hinge_analysis.mat'])
        plot_dir = [output_dir filesep 'Hinge_Plots'];
        for i = 1:height(hinge)
            plot(hinge.rotation_TH{i},hinge.moment_TH{i});
            ylabel('Hinge Moment (k-in)')
            xlabel('Hinge Rotation (rads)')
            plot_name = ['Hinge ' num2str(i) ' Response'];
            fn_format_and_save_plot( plot_dir, plot_name, 2 )
        end
    end
%% Plot PM Diagrams for each element
%     for i = 1:length(element.id)
%         if strcmp(element.type{i},'column')
%             ele = element(i,:);
%             ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
%             ele_PM = element_PM.(['ele_' num2str(element.id(i))]);
% 
%             hold on
%             plot(ele_PM.vector_M/1000,ele_PM.vector_P/1000,'k','LineWidth',2)
%             plot(abs(ele_TH.M_TH_1)/1000,ele_TH.P_TH_1/1000,'b','LineWidth',0.75)
%             ylabel('Axial (k)')
%             xlabel('Moment (k-in)')
%             plot_dir = [output_dir filesep 'PM_plots'];
%             plot_name = ['ele_' num2str(element.id(i))];
%             fn_format_and_save_plot( plot_dir, plot_name, 2 )
%         end
% 
%     end
end

