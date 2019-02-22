function [ ] = fn_plot_hinge_response( plot_dir, hinge, element, ele_prop_table, node, stories2plot, eq_analysis_timespace )
% Description: Fn to plot hinge backbone along with hinge response.

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import plotting_tools.fn_format_and_save_plot
import plotting_tools.fn_plot_backbone
import asce_41.fn_define_backbone_rot


% Define plot directory
hinge_plot_dir = [plot_dir filesep 'Hinge_Plots'];

%% Begin Method
for i = 1:height(hinge)
    % Grab Element Properties
    ele = element(element.id == hinge.element_id(i),:);
    if ele.story <= stories2plot
        ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);

        % Define plot directory
        hinge_plot_dir = [plot_dir filesep 'Hinge_Plots' filesep 'Story - ' num2str(ele.story)];

        % Set critical mode variable
        if strcmp(hinge.direction{i},'oop')
            crit_mode = ele.critical_mode_oop;
        else
            crit_mode = ele.critical_mode;
        end

        if strcmp(hinge.type(i),'rotational')
            if strcmp(ele.type,'column')
                hinge_name = ['Hinge_y_', num2str(node.y(node.id == hinge.node_1(i))),'_',hinge.direction{i}];
            elseif strcmp(ele.type,'beam')
                hinge_name = ['Hinge_x_', num2str(node.x(node.id == hinge.node_1(i))),'_',hinge.direction{i}];
            end
            plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - ' hinge_name ' - Rotation Response'];


            fn_plot_backbone( ele, ele_props, hinge_plot_dir, plot_name, 2, hinge.deformation_TH{i}, hinge.force_TH{i}, crit_mode, hinge.direction{i})

            % Plot Hinge Rotation Time History
            [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.Mn_pos, ele.Mn_neg, ele.Mp_pos, ele.Mp_neg, ele.length, ele_props.e, ele_props.iz, ele.a_hinge, ele.b_hinge, ele.c_hinge, 10, 0.1 );
            hold on
            pos_yeild = rot_vec_pos(1);
            neg_yeild = rot_vec_neg(1);
            pos_a = ele.a_hinge + pos_yeild;
            neg_a = ele.a_hinge + neg_yeild;
            pos_b = ele.b_hinge + pos_yeild;
            neg_b = ele.b_hinge + neg_yeild;
            yeild_time_point = min(eq_analysis_timespace(hinge.deformation_TH{i}>=pos_yeild | hinge.deformation_TH{i}<=-neg_yeild));
            a_time_point = min(eq_analysis_timespace(hinge.deformation_TH{i}>=pos_a | hinge.deformation_TH{i}<=-neg_a));
            b_time_point = min(eq_analysis_timespace(hinge.deformation_TH{i}>=pos_b | hinge.deformation_TH{i}<=-neg_b));
            max_deform = 1.25*max(abs(hinge.deformation_TH{i}));
            if ~isempty(yeild_time_point)
                hist_plot = plot([yeild_time_point,yeild_time_point],[-max_deform,max_deform],'--k','LineWidth',1.25,'DisplayName','Yield');
            end
            if ~isempty(a_time_point)
                hist_plot = plot([a_time_point,a_time_point],[-max_deform,max_deform],'--m','LineWidth',1.25,'DisplayName','Ultimate Capacity');
            end
            if ~isempty(b_time_point)
                hist_plot = plot([b_time_point,b_time_point],[-max_deform,max_deform],'--r','LineWidth',1.25,'DisplayName','Failure');
            end
            hist_plot = plot(eq_analysis_timespace,hinge.deformation_TH{i},'b','LineWidth',1,'DisplayName','Analysis');
            ylabel('Hinge Rotation (rads)')
            xlabel('Time (s)')
            xlim([0,15])
            ylim([-max_deform,max_deform])
            plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - ' hinge_name ' - Time History'];
            fn_format_and_save_plot( hinge_plot_dir, plot_name, 1 )

        elseif strcmp(hinge.type(i),'shear')
            plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - Shear Response'];
            fn_plot_backbone( ele, ele_props, hinge_plot_dir, plot_name, 2, hinge.deformation_TH{i}, hinge.force_TH{i}, crit_mode, hinge.direction{i})
        end
    end
end
end

