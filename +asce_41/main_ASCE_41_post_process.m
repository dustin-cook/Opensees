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
load([read_dir filesep 'hinge_analysis.mat'])
load([read_dir filesep 'node_analysis.mat'])

% mf_joint_table = readtable([analysis.out_dir filesep 'model_data' filesep 'mf_joint.csv'],'ReadVariableNames',true);
% joint_table = readtable([analysis.out_dir filesep 'model_data' filesep 'joint.csv'],'ReadVariableNames',true);

%% Calculate Element Properties and Modify Analysis Results based on ASCE 41-17
% Basic building or analysis properties
[ model, element, torsion, node ] = fn_basic_analysis_properties( model, story, element, node );

if analysis.asce_41_post_process
    % Procedure Specific Analysis
    if strcmp(analysis.proceedure,'NDP') % Nonlinear Dynamic Proceedure
        if analysis.nonlinear ~= 0 && analysis.type == 1 % Only for last nonlinear dynamic run
            OS_demands = element;
            OS_joint = joint;
%             if length(analysis.type_list) > 1
                load([write_dir filesep 'element_analysis.mat'])
                load([write_dir filesep 'joint_analysis.mat'])
                load([write_dir filesep 'hinge_analysis.mat'])
%             else
%                 load([read_dir filesep 'element_analysis.mat'])
%                 load([read_dir filesep 'joint_analysis.mat'])
%                 load([read_dir filesep 'hinge_analysis.mat'])
%             end
            % Get demands from recent runs
            element.Pmax = OS_demands.Pmax;
            element.Pmin = OS_demands.Pmin;
            element.Vmax_1 = OS_demands.Vmax_1;
            element.Vmax_2 = OS_demands.Vmax_2;
            element.Vmax_oop_1 = OS_demands.Vmax_oop_1;
            element.Vmax_oop_2 = OS_demands.Vmax_oop_2;
            element.Mmax_1 = OS_demands.Mmax_1;
            element.Mmax_2 = OS_demands.Mmax_2;
            [ OS_joint ] = fn_joint_capacity( OS_joint, element, ele_prop_table );
            joint.Pmax = OS_joint.Pmax;
            joint.Vmax = OS_joint.Vmax;
            % calculate hinge demand to capacity ratios
            [ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge, read_dir, node );
            
            % Filter Accelerations
            [ story ] = fn_accel_filter(node, story, analysis.filter_freq_range, read_dir, write_dir);
            
        else % Pushover runs
            [ element, joint ] = main_element_capacity( story, ele_prop_table, element, analysis, joint, read_dir, write_dir );
            [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );
        end
    elseif strcmp(analysis.proceedure,'LDP') % Linear Dynamic Proceedure and Test Proceedure (and all others defined so be careful)
        fn_torsional_amplification( story, element ) % Torsion check (will throw errors if triggered, need to update)
        [ model, element, joint ] = fn_linear_capacity_and_c_factors( model, story, ele_prop_table, element, analysis, joint, read_dir, write_dir );
        [ element ] = main_m_factors( ele_prop_table, element );
        [ element, ~ ] = fn_calc_dcr( element, element_TH, 'cp' );
    else % Test Cases
        [ element, joint ] = main_element_capacity( story, ele_prop_table, element, analysis, joint, read_dir, write_dir );
        [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );
    end
end

%% Save Data
save([write_dir filesep 'model_analysis.mat'],'model')
save([write_dir filesep 'story_analysis.mat'],'story')
save([write_dir filesep 'element_analysis.mat'],'element')
save([write_dir filesep 'hinge_analysis.mat'],'hinge')
save([write_dir filesep 'joint_analysis.mat'],'joint')
save([write_dir filesep 'node_analysis.mat'],'node')
% Save load case info
if ~strcmp(analysis.case,'NA')
    write_dir = [analysis.out_dir filesep analysis.case];
    fn_make_directory( write_dir )
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
end

% Save capacities to compare iterations
capacity = element.capacity;
end

