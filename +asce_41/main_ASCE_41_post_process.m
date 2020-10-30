function [ capacity, torsion ] = main_ASCE_41_post_process( analysis, ele_prop_table )
% Description: Main script that post process an ASCE 41 analysis

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import asce_41.*

% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'opensees_data'];
write_dir = [analysis.out_dir filesep 'asce_41_data'];
if ~exist(write_dir,'dir')
    fn_make_directory( write_dir )
end

% Load Analysis Data
load([read_dir filesep 'model_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'joint_analysis.mat'])
load([read_dir filesep 'node_analysis.mat'])

%% Calculate Element Properties and Modify Analysis Results based on ASCE 41-17
% Basic building or analysis properties
[ model, element, torsion, node ] = fn_basic_analysis_properties( model, story, element, node );

% Filter Accelerations
if analysis.type == 1 && analysis.filter_accel == 1
    [ story ] = fn_accel_filter(node, story, analysis.filter_freq_range, read_dir, write_dir);
end

if analysis.asce_41_post_process
    % Procedure Specific Analysis
    if strcmp(analysis.proceedure,'NDP') && analysis.type == 1 % Nonlinear Dynamic Proceedure
        % Merge demands from analysis into capacities from pushover
        [element] = fn_merge_demands(element, write_dir, model);
        [joint] = fn_merge_joint_demands(joint, element, ele_prop_table, write_dir, read_dir);
        
        % calculate hinge demand to capacity ratios
        load([read_dir filesep 'hinge_analysis.mat'])
        [ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge, read_dir, node );
        [ joint ] = fn_accept_joint( joint, analysis.joint_explicit, read_dir);
            
    elseif strcmp(analysis.proceedure,'LDP') && analysis.type == 1 % Linear Dynamic Proceedure
        % Merge demands from analysis into capacities from pushover
        [element] = fn_merge_demands(element, write_dir, model);
        [joint] = fn_merge_joint_demands(joint, element, ele_prop_table, write_dir);
      
    else % Pushovers and all other cases
        [ element, joint ] = main_element_capacity( story, ele_prop_table, element, analysis, joint, read_dir, write_dir );
        [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );
    end
end

%% Save Data
save([write_dir filesep 'model_analysis.mat'],'model')
save([write_dir filesep 'story_analysis.mat'],'story')
save([write_dir filesep 'element_analysis.mat'],'element')
save([write_dir filesep 'joint_analysis.mat'],'joint')
save([write_dir filesep 'node_analysis.mat'],'node')
if exist('hinge','var')
    save([write_dir filesep 'hinge_analysis.mat'],'hinge')
end

% Save load case info
if ~strcmp(analysis.case,'NA')
    write_dir = [analysis.out_dir filesep analysis.case];
    fn_make_directory( write_dir )
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
    if exist('hinge','var')
        save([write_dir filesep 'hinge_analysis.mat'],'hinge')
    end
end

% % Save capacities to compare iterations
% if analysis.asce_41_post_process
%     capacity = element.capacity;
% else
%     capacity = [];
% end
end

function [element] = fn_merge_demands(element, write_dir, model)

% Save demand tables from most recent opensees run
OS_demands = element;

% Load Post process capacity data
load([write_dir filesep 'element_analysis.mat'])

% Merge demands into capacity tables
element.Pmax = OS_demands.Pmax;
element.Pmin = OS_demands.Pmin;
element.Vmax_1 = OS_demands.Vmax_1;
element.Vmax_2 = OS_demands.Vmax_2;
if strcmp(model.dimension,'3D')
    element.Vmax_oop_1 = OS_demands.Vmax_oop_1;
    element.Vmax_oop_2 = OS_demands.Vmax_oop_2;
end
element.Mmax_1 = OS_demands.Mmax_1;
element.Mmax_2 = OS_demands.Mmax_2;

% Keep node info from latest run (for linear sake)
element.node_1 = OS_demands.node_1;
element.node_2 = OS_demands.node_2;

end

function [joint] = fn_merge_joint_demands(joint, element, ele_prop_table, write_dir, read_dir)
%% Joints
% import packages
import asce_41.fn_joint_capacity

% Load Post process capacity data
load([write_dir filesep 'joint_analysis.mat'])

% Save demand tables from most recent opensees run
OS_joint = joint;

% Determine Joint demands based on element capacities
for i = 1:height(OS_joint)
    OS_jnt = OS_joint(i,:);
    TH_file = [read_dir filesep 'joint_TH_' num2str(OS_jnt.id) '.mat'];
    if ~exist(TH_file,'file')
        jnt_TH = [];
    else
        load(TH_file)
    end
    [ OS_jnt ] = fn_joint_capacity( OS_jnt, element, ele_prop_table, jnt_TH );

    % Merge joint demands into capacity table
    joint.Pmax(i) = OS_jnt.Pmax;
    joint.Vmax(i) = OS_jnt.Vmax;
end


end
