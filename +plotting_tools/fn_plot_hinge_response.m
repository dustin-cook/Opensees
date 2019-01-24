function [ ] = fn_plot_hinge_response( plot_dir, hinge, element, ele_prop_table, node )
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

% Define plot directory
hinge_plot_dir = [plot_dir filesep 'Hinge_Plots'];

%% Begin Method
for i = 1:height(hinge)
    % Grab Element Properties
    ele = element(element.id == hinge.element_id(i),:);
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

%         % Plot Hinge Rotation Time History
%         hold on
%         yeild_point = theta_yeild - (10/11)*theta_yeild;
%         b_point = ele.b_hinge + yeild_point;
%         hist_plot = plot([0,15],[yeild_point,yeild_point],'--','color',[0.5,0.5,0.5],'LineWidth',1.25,'DisplayName','yield');
%         set(get(get(hist_plot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%         hist_plot = plot([0,15],[b_point,b_point],'--k','LineWidth',1.25,'DisplayName','b');
%         set(get(get(hist_plot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%         hist_plot = plot([0,15],[-yeild_point,-yeild_point],'--','color',[0.5,0.5,0.5],'LineWidth',1.25,'DisplayName','yield');
%         hist_plot = plot([0,15],[-b_point,-b_point],'--k','LineWidth',1.25,'DisplayName','b');
%         hist_plot = plot(eq_analysis_timespace,hinge.rotation_TH{i},'b','LineWidth',1,'DisplayName','Analysis');
%         ylabel('Hinge Rotation (rads)')
%         xlabel('Time (s)')
%         xlim([0,15])
%         ylim([-1.5*b_point,1.5*b_point])
%         plot_name = ['element_' num2str(hinge.element_id(i)) ' - ' hinge_name ' Rotation Time History'];
%         fn_format_and_save_plot( hinge_plot_dir, plot_name, 2 )

    elseif strcmp(hinge.type(i),'shear')
        plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - Shear Response'];
        fn_plot_backbone( ele, ele_props, hinge_plot_dir, plot_name, 2, hinge.deformation_TH{i}, hinge.force_TH{i}, crit_mode, hinge.direction{i})
    end
end
end

