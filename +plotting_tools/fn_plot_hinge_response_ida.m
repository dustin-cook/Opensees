function [ ] = fn_plot_hinge_response_ida( plot_dir, read_dir, hinge, element, ele_prop_table, node, analysis )
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
% hinge_plot_dir = [plot_dir filesep 'Hinge_Plots'];

%% Plot Element Hinges
for i = 1:height(hinge)
    % Grab Element Properties
    ele = element(element.id == hinge.element_id(i),:);
    ele_side = num2str(hinge.ele_side(i));
    if ele.story <= analysis.hinge_stories_2_plot
        if exist([read_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'file')
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
                elseif strcmp(ele.type,'wall')
                    hinge_name = ['Hinge_wall_', num2str(node.x(node.id == hinge.node_1(i))),'_',hinge.direction{i}];
                end
                plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - ' hinge_name ' - Rotation Response'];
                fn_plot_backbone( ele, ele_side, ele, read_dir, hinge_plot_dir, plot_name, 2, hin_TH.deformation_TH, hin_TH.force_TH, crit_mode, hinge.direction{i})

            elseif strcmp(hinge.type(i),'shear')
                plot_name = [ele.type{1} '_' num2str(hinge.element_id(i)) ' - Shear Response'];
                fn_plot_backbone( ele, ele_side, ele_props, read_dir, hinge_plot_dir, plot_name, 2, hin_TH.deformation_TH, hin_TH.force_TH, crit_mode, hinge.direction{i})
            end
        end
    end
end


end

