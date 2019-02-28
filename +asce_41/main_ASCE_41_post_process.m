function [ capacity ] = main_ASCE_41_post_process( analysis, ele_prop_table )
% Description: Main script that post process an ASCE 41 analysis

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import asce_41.main_element_capacity
import asce_41.main_m_factors
import asce_41.fn_linear_capacity_and_c_factors
import asce_41.main_hinge_properties
import asce_41.fn_accept_hinge
import asce_41.fn_torsional_amplification
import asce_41.fn_basic_analysis_properties
import asce_41.fn_calc_dcr

% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'opensees_data'];
write_dir = [analysis.out_dir filesep 'asce_41_data'];
fn_make_directory( write_dir )

% Load Analysis Data
load([read_dir filesep 'model_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'joint_analysis.mat'])
load([read_dir filesep 'element_TH.mat'])
load([read_dir filesep 'hinge_analysis.mat'])

% mf_joint_table = readtable([analysis.out_dir filesep 'model_data' filesep 'mf_joint.csv'],'ReadVariableNames',true);
% joint_table = readtable([analysis.out_dir filesep 'model_data' filesep 'joint.csv'],'ReadVariableNames',true);

%% Calculate Element Properties and Modify Analysis Results based on ASCE 41-17
% Basic building or analysis properties
[ model, element ] = fn_basic_analysis_properties( model, story, element );

if analysis.asce_41_post_process
    % Torsion Check and Amplification
    fn_torsional_amplification( story, element ) % NEED TO UPDATE

    % Procedure Specific Analysis
    if strcmp(analysis.proceedure,'NDP') % Nonlinear Dynamic Proceedure
        [ element, element_TH, element_PM, joint ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis, joint  );
        [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );
        if analysis.nonlinear ~= 0 % Only for nonlinear runs
            [ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge );
        end
    elseif strcmp(analysis.proceedure,'LDP') % Linear Dynamic Proceedure and Test Proceedure (and all others defined so be careful)
        [ model, element, element_TH, element_PM, joint ] = fn_linear_capacity_and_c_factors( model, story, ele_prop_table, element, element_TH, analysis, joint );
        [ element ] = main_m_factors( ele_prop_table, element );
        [ element, ~ ] = fn_calc_dcr( element, element_TH, 'cp' );
    else % Test Cases
        [ element, element_TH, element_PM, joint ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis, joint  );
        [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint );
    end
    
    %% Save Data
    save([write_dir filesep 'element_PM.mat'],'element_PM')
end

%% Save Data
save([write_dir filesep 'model_analysis.mat'],'model')
save([write_dir filesep 'story_analysis.mat'],'story')
save([write_dir filesep 'element_analysis.mat'],'element')
save([write_dir filesep 'element_TH.mat'],'element_TH')
save([write_dir filesep 'hinge_analysis.mat'],'hinge')
save([write_dir filesep 'joint_analysis.mat'],'joint')

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

