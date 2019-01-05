function [ ] = main_check_analysis( analysis, ele_prop_table )
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
write_dir = [analysis.out_dir filesep 'validation_plots'];
fn_make_directory( write_dir )

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])

%% Plot Hinge Convergence
if analysis.plot_hinges && strcmp(analysis.proceedure,'NDP')
    for i = 1:height(element)
        ele = element(i,:);
        ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        plot_name = ['element_' num2str(ele.id)];
        fn_plot_backbone( ele, ele_props, write_dir, plot_name, 1)
    end
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

%% Calculate Beam Column Strength Ratios
% for i =1:length(joint.id)
%    beam1 = max([element.Mn_pos(element.node_2 == joint.x_neg(i)),element.Mn_neg(element.node_2 == joint.x_neg(i))]); % Maximum of beam pos and neg nominal bending strength
%    beam2 = max([element.Mn_pos(element.node_1 == joint.x_pos(i)),element.Mn_neg(element.node_1 == joint.x_pos(i))]); 
%    column1 = min([element.Mn_pos(element.node_2 == joint.y_neg(i)),element.Mn_neg(element.node_2 == joint.y_neg(i))]); % Minimum of column pos and negative nominal moment strength
%    column2 = min([element.Mn_pos(element.node_1 == joint.y_pos(i)),element.Mn_neg(element.node_1 == joint.y_pos(i))]); 
%    joint.beam_strength(i) = sum([beam1,beam2]);
%    joint.column_strength(i) = sum([column1,column2]);
%    joint.col_bm_ratio(i) = joint.column_strength(i)/joint.beam_strength(i);
% end

end

