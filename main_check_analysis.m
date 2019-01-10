function [ ] = main_check_analysis( analysis, ele_prop_table, capacity, step )
% Description: Checks proceedures and analyses are working as expected and
% creates visuals. Check come from both a general modeling perpective as
% well as specific checks perscribed in ASCE 41-17.

% Created By: Dustin Cook
% Date Created: 1/4/2019

% Inputs: 

% Outputs: 

% Assumptions:

%% Initial Setup
% Import Packages
import plotting_tools.fn_plot_backbone

% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'asce_41_data'];
write_dir = [analysis.out_dir filesep 'validation_plots' filesep 'hinge_plots'];

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])

%% Plot Hinge Convergence
if analysis.element_plots && strcmp(analysis.proceedure,'NDP')
    for i = 1:height(element)

        ele = element(i,:);
        ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        plot_name = ['element_' num2str(ele.id)];
        
        prev_fig_file = [write_dir filesep plot_name '.fig'];
        if exist(prev_fig_file,'file')
            openfig(prev_fig_file);
            hold on
        end
        
        line_color = [1,1,1] - step/length(analysis.type_list);
        fn_plot_backbone( ele, ele_props, write_dir, plot_name, 1, 0, 0, line_color)
    end
end

%% Check Capacity Convergence
if step == length(analysis.type_list)
    capacity
    plot([1,2,3,4],capacity)
end

%% Vertical Ground Motion Convergence


%% Torsion
% If the displacement multiplier ? caused by actual plus accidental 
% torsion at any level exceeds 1.5, two-dimensional models shall not be 
% permitted and three-dimensional models that account for the spatial 
% distribution of mass and stiffness shall be used.

% Increased forces and displacements caused by accidental torsion need not 
% be considered if either of the following conditions apply: (a) the 
% accidental torsional moment is less than 25% of the actual torsional 
% moment, or (b) the ratio of the displacement multiplier ? caused by the 
% actual plus accidental torsion and the displacement multiplier caused by 
% actual torsion is less than 1.1 at every floor.

end
