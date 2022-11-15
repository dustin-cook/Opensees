function [ ] = fn_plot_hinge_response_ida_fiber( plot_dir, read_dir, element, ele_prop_table, analysis )
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
for i = 1:height(element)
    % Grab Element Properties
    ele = element(i,:);
    ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    if ele.story <= analysis.hinge_stories_2_plot
        if exist([read_dir filesep 'element_TH_' num2str(element.id(i)) '.mat'],'file')
            load([read_dir filesep 'element_TH_' num2str(element.id(i)) '.mat'])

            % Define plot directory
            hinge_plot_dir = [plot_dir filesep 'Hinge_Plots' filesep 'Story - ' num2str(ele.story)];

            % Plot backbone
            plot_name = [ele.type{1} '_' num2str(ele.id) ' - Rotation Response'];
%             if strcmp(ele.type,'column')
%                 M_sign = -1;
%             else
                M_sign = 1;
%             end
            fn_plot_backbone( ele, [], ele_prop, read_dir, hinge_plot_dir, plot_name, 2, ele_TH.rot, M_sign*ele_TH.M_TH_1, [], [])
        end
    end
end


end

