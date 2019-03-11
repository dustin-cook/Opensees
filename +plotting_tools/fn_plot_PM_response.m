function [ ] = fn_plot_PM_response( plot_dir, read_dir, element, stories_2_plot )
% Description: Fn to plot PM interaction diagrams along with column
% response.

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Begin Method
for i = 1:length(element.id)
    ele = element(i,:);
    if ele.story <= stories_2_plot
        load([read_dir filesep 'element_TH_' num2str(ele.id) '.mat'])
        if strcmp(element.type{i},'column')
            load([read_dir filesep 'element_PM_' num2str(ele.id) '.mat'])

            hold on
            plot(ele_PM.vector_M_1/1000,ele_PM.vector_P_1/1000,'k','LineWidth',2) % Only plotting at the bottom hinge for now
            plot(abs(ele_TH.M_TH_1)/1000,ele_TH.P_TH_1/1000,'b','LineWidth',0.75)
            ylabel('Axial (k)')
            xlabel('Moment (k-in)')
            pm_plot_dir = [plot_dir filesep 'PM Plots' filesep 'Story - ' num2str(ele.story)];
            plot_name = ['column_' num2str(element.id(i))];
            fn_format_and_save_plot( pm_plot_dir, plot_name, 2 )
        end
    end
end

end

