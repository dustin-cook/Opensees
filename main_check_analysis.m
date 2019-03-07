function [ ] = main_check_analysis( analysis, ele_prop_table, capacity, torsion, step )
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
import plotting_tools.fn_format_and_save_plot

% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'asce_41_data'];
write_dir = [analysis.out_dir filesep 'validation_plots'];

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'joint_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])

% Define Analysis Checklist File
file_name = [analysis.out_dir filesep 'analysis_checklist.txt'];
if step == 1
    fileID = fopen(file_name,'w');
else
    fileID = fopen(file_name,'a');
end

%% Plot Hinge Convergence
if analysis.element_plots && strcmp(analysis.proceedure,'NDP') && analysis.type == 2
    for i = 1:height(element)
        ele = element(i,:);
        if ele.story <= analysis.hinge_stories_2_plot
            ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
            plot_dir = [write_dir filesep 'hinge_plots' filesep 'Story - ' num2str(ele.story)];
            plot_name = [ele.type{1} '_' num2str(ele.id)];

            prev_fig_file = [plot_dir filesep plot_name '.fig'];
            if exist(prev_fig_file,'file')
                openfig(prev_fig_file);
                hold on
            end

            line_color = [1,1,1] - step/(length(analysis.type_list)-1);
            fn_plot_backbone( ele, '1', ele_props, plot_dir, plot_name, 1, 0, 0, ele.critical_mode_1, 'primary', line_color)
        end
    end
end

%% Check Capacity Convergence
if strcmp(analysis.case,'backbones') && strcmp(analysis.proceedure,'NDP')
    relative_capacity = capacity(:,2:end) ./ capacity(:,1:(end-1));
    iter_num = 2:(length(relative_capacity(1,:))+1);
    plot(iter_num,relative_capacity)
    ylim([0,2])
    ylabel('Percent Change in Capacity')
    xlabel('Iteration')
    plot_name = 'Capacity Convergence';
    fn_format_and_save_plot( write_dir, plot_name, 2 )
end

%% Joint Shear Demand is Less than Capacity
if strcmp(analysis.proceedure,'NDP') && analysis.type == 1 && analysis.nonlinear == 1
    joint_shear_ratio = joint.Vmax ./ joint.Vj;
    fprintf(fileID,'Max Joint Shear Ratio =  %f\n',max(joint_shear_ratio));
end

%% Walls do not yeild
if strcmp(analysis.proceedure,'NDP') && analysis.type == 1 && analysis.nonlinear == 1
    walls = element(strcmp(element.type,'wall'),:);
    wall_shear_ratio = walls.Vmax ./ (walls.f_hinge_1 .* walls.Vn_1);
    fprintf(fileID,'Max Wall Shear Ratio =  %f of first yield (f-value)\n',max(wall_shear_ratio));
end

%% Vertical Ground Motion Convergence


%% Torsion
% If the displacement multiplier eta caused by actual plus accidental 
% torsion at any level exceeds 1.5, two-dimensional models shall not be 
% permitted and three-dimensional models that account for the spatial 
% distribution of mass and stiffness shall be used.

% Increased forces and displacements caused by accidental torsion need not 
% be considered if either of the following conditions apply: (a) the 
% accidental torsional moment is less than 25% of the actual torsional 
% moment, or (b) the ratio of the displacement multiplier eta caused by the 
% actual plus accidental torsion and the displacement multiplier caused by 
% actual torsion is less than 1.1 at every floor.

if strcmp(analysis.case,'torsion_check') && analysis.type == 1
    % Case b check
    eta_ratio_x = torsion{step}.x./torsion{step-1}.x;
    eta_ratio_z = torsion{step}.z./torsion{step-1}.z;
    max_eta_ratio = max([eta_ratio_x; eta_ratio_z]);
    fprintf(fileID,'Max Accidental Torsion Effect = %f of Actual Torsion (f-value)\n',max_eta_ratio);
    
    % 2D applicability check
    max_eta = max([torsion{step}.x;torsion{step}.z]);
    if max_eta > 1.5 
        fprintf(fileID,'2D models are not allowed: Max Accidental and Actual Torsion Factor = %f \n',max_eta);
    end
end



% Close File
fclose(fileID);

end

