function [ ] = fn_plot_PM_response( plot_dir, element, element_TH, element_PM )
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
    ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
    if strcmp(element.type{i},'column')
        ele_PM = element_PM.(['ele_' num2str(element.id(i))]);

        hold on
        plot(ele_PM.vector_M/1000,ele_PM.vector_P/1000,'k','LineWidth',2)
        plot(abs(ele_TH.M_TH_1)/1000,ele_TH.P_TH_1/1000,'b','LineWidth',0.75)
        ylabel('Axial (k)')
        xlabel('Moment (k-in)')
        pm_plot_dir = [plot_dir filesep 'PM Plots' filesep 'Story - ' num2str(ele.story)];
        plot_name = ['column_' num2str(element.id(i))];
        fn_format_and_save_plot( pm_plot_dir, plot_name, 2 )
    end
end

end

