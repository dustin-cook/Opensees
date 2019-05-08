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
os_dir = [analysis.out_dir filesep 'opensees_data'];
read_dir = [analysis.out_dir filesep 'asce_41_data'];
write_dir = [analysis.out_dir filesep 'validation_plots'];

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'joint_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'hinge_analysis.mat'])

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
if analysis.type == 1 && strcmp(analysis.proceedure,'NDP') % && strcmp(analysis.case,'backbones')
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
    if ~isempty(walls)
        wall_shear_ratio = walls.Vmax_1 ./ (walls.f_hinge_1 .* walls.Vn_1);
        fprintf(fileID,'Max Wall Shear Ratio =  %f of first yield (f-value)\n',max(wall_shear_ratio));
    end
end

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
    fprintf(fileID,'Max Accidental Torsion Effect = %f of Actual Torsion \n',max_eta_ratio);
    
    % 2D applicability check
    max_eta = max([torsion{step}.x;torsion{step}.z]);
    if max_eta > 1.5 
        fprintf(fileID,'2D models are not allowed: Max Accidental and Actual Torsion Factor = %f \n',max_eta);
    end
end

%% Equilibrium of hinge forces with element forces
if analysis.element_plots && analysis.nonlinear == 1 && analysis.type == 1
    for i = 1:height(hinge)
        hin = hinge(i,:);
        ele = element(element.id == hin.element_id,:);
        if ele.story <= analysis.hinge_stories_2_plot
            % Load time history data
            load([os_dir filesep 'element_TH_' num2str(ele.id) '.mat'])
            load([os_dir filesep 'hinge_TH_' num2str(hin.id) '.mat'])
            
            
            % Plot Element v Hinge Shear
            if strcmp(ele.type,'wall') && ~strcmp(hin.direction,'oop')
                hold on
                plot(ele_TH.(['V_TH_' num2str(hin.ele_side)]),'DisplayName','Element Force Recorder')
                plot(hin_TH.shear_TH,'--k','DisplayName','Element Hinge Recorder')
                ylabel('Shear Force (lbs)')
                plot_name = [ele.type{1} ' ' num2str(ele.id) ' side ' num2str(hin.ele_side) ' - ' hin.direction{1} ' shear']; 
                fn_format_and_save_plot( write_dir, plot_name, 1 )
            else    
                % Plot Element v Hinge Moments
                hold on
                if strcmp(hin.direction,'oop')
                    plot(ele_TH.(['M_TH_oop_' num2str(hin.ele_side)]),'DisplayName','Element Force Recorder')
                else
                    plot(ele_TH.(['M_TH_' num2str(hin.ele_side)]),'DisplayName','Element Force Recorder')
                end
                plot(hin_TH.moment_TH,'--k','DisplayName','Element Hinge Recorder')
                ylabel('Moment Demand (lbs-in)')
                plot_name = [ele.type{1} ' ' num2str(ele.id) ' side ' num2str(hin.ele_side) ' - ' hin.direction{1} ' moment']; 
                fn_format_and_save_plot( write_dir, plot_name, 1 )
            end
            
        end
    end
end

%% Vertical Ground Motion Convergence

% Close File
fclose(fileID);

end

