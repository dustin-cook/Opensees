function [ ] = fn_plot_hinge_response( plot_dir, read_dir, hinge, element, ele_prop_table, node, stories2plot, eq_analysis_timespace )
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
    ele_side = num2str(hinge.ele_side(i));
    if ele.story <= stories2plot
        ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        load([read_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'])
        
        % Define plot directory
        hinge_plot_dir = [plot_dir filesep 'Hinge_Plots' filesep 'Story - ' num2str(ele.story)];

        % Set critical mode variable
        if strcmp(hinge.direction{i},'oop')
            crit_mode = ele.(['critical_mode_oop_' ele_side]);
        else
            crit_mode = ele.(['critical_mode_' ele_side]);
        end

        if strcmp(hinge.type(i),'rotational')
            if strcmp(ele.type,'column')
                hinge_name = ['Hinge_y_', num2str(node.y(node.id == hinge.node_1(i))),'_',hinge.direction{i}];
            elseif strcmp(ele.type,'beam')
                hinge_name = ['Hinge_x_', num2str(node.x(node.id == hinge.node_1(i))),'_',hinge.direction{i}];
            end
            plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - ' hinge_name ' - Rotation Response'];


            fn_plot_backbone( ele, ele_side, ele_props, read_dir, hinge_plot_dir, plot_name, 2, hin_TH.deformation_TH, hin_TH.force_TH, crit_mode, hinge.direction{i})

            % Plot Hinge Rotation Time History
%             if exist('eq_analysis_timespace','var')
%                 [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.(['Mn_pos_' ele_side]), ele.(['Mn_neg_' ele_side]), ele.(['Mp_pos_' ele_side]), ele.(['Mp_neg_' ele_side]), ele.length, ele_props.e, ele_props.iz, ele.(['a_hinge_' ele_side]), ele.(['b_hinge_' ele_side]), ele.(['c_hinge_' ele_side]), 10, 0.1, crit_mode );
%                 hold on
%                 pos_yeild = rot_vec_pos(1);
%                 neg_yeild = rot_vec_neg(1);
%                 pos_a = ele.(['a_hinge_' ele_side]) + pos_yeild;
%                 neg_a = ele.(['a_hinge_' ele_side]) + neg_yeild;
%                 pos_b = ele.(['b_hinge_' ele_side]) + pos_yeild;
%                 neg_b = ele.(['b_hinge_' ele_side]) + neg_yeild;
%                 yeild_time_point = min(eq_analysis_timespace(hin_TH.deformation_TH>=pos_yeild | hin_TH.deformation_TH<=-neg_yeild));
%                 a_time_point = min(eq_analysis_timespace(hin_TH.deformation_TH>=pos_a | hin_TH.deformation_TH<=-neg_a));
%                 b_time_point = min(eq_analysis_timespace(hin_TH.deformation_TH>=pos_b | hin_TH.deformation_TH<=-neg_b));
%                 max_deform = 1.25*max(abs(hin_TH.deformation_TH));
%                 if ~isempty(yeild_time_point)
%                     hist_plot = plot([yeild_time_point,yeild_time_point],[-max_deform,max_deform],'--k','LineWidth',1.25,'DisplayName','Yield');
%                 end
%                 if ~isempty(a_time_point)
%                     hist_plot = plot([a_time_point,a_time_point],[-max_deform,max_deform],'--m','LineWidth',1.25,'DisplayName','Ultimate Capacity');
%                 end
%                 if ~isempty(b_time_point)
%                     hist_plot = plot([b_time_point,b_time_point],[-max_deform,max_deform],'--r','LineWidth',1.25,'DisplayName','Failure');
%                 end
%                 hist_plot = plot(eq_analysis_timespace,hin_TH.deformation_TH,'b','LineWidth',1,'DisplayName','Analysis');
%                 ylabel('Hinge Rotation (rads)')
%                 xlabel('Time (s)')
%                 xlim([0,15])
%                 ylim([-max_deform,max_deform])
%                 plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - ' hinge_name ' - Time History'];
%                 fn_format_and_save_plot( hinge_plot_dir, plot_name, 1 )
%             end

        elseif strcmp(hinge.type(i),'shear')
            plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - Shear Response'];
            fn_plot_backbone( ele, ele_side, ele_props, hinge_plot_dir, plot_name, 2, hin_TH.deformation_TH, hin_TH.force_TH, crit_mode, hinge.direction{i})
        end
    end
end
end

