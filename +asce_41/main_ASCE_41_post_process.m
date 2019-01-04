function [ ] = main_ASCE_41_post_process( analysis )
% Description: Main script that post process an ASCE 41 analysis

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:

%% Import Packages
import asce_41.main_element_capacity
import asce_41.main_m_factors
import asce_41.fn_linear_capacity_and_c_factors
import asce_41.main_hinge_properties
import asce_41.fn_accept_hinge
import asce_41.fn_torsional_amplification
import asce_41.fn_basic_analysis_properties
import asce_41.fn_calc_dcr

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
m_table.col = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.beam = readtable(['+asce_41' filesep 'linear_beam_m.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'model_analysis.mat'])
load([output_dir filesep 'story_analysis.mat'])
load([output_dir filesep 'element_analysis.mat'])
load([output_dir filesep 'element_TH.mat'])
load([output_dir filesep 'hinge_analysis.mat'])

%% Calculate Element Properties and Modify Analysis Results based on ASCE 41-17
% Basic building or analysis properties
[ model, element ] = fn_basic_analysis_properties( model, story, element );

% Torsion Check and Amplification
fn_torsional_amplification( story, element ) % NEED TO UPDATE

% Procedure Specific Analysis
if strcmp(analysis.proceedure,'LDP') % Linear Dynamic Proceedure
    [ model, element, element_TH, element_PM ] = fn_linear_capacity_and_c_factors( model, story, ele_prop_table, element, element_TH, analysis );
    [ element ] = main_m_factors( ele_prop_table, element, m_table );
    [ element, ~ ] = fn_calc_dcr( element, element_TH, 'cp' );
elseif strcmp(analysis.proceedure,'NDP') % Nonlinear Dynamic Proceedure
    [ element, element_TH, element_PM ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis );
    if analysis.nonlinear == 0 % First Linear Run to get hinge parameters
        [ element ] = main_hinge_properties( ele_prop_table, element, 0 );
    else % Rest of the nonlinear runs
        [ element ] = main_hinge_properties( ele_prop_table, element, analysis.plot_hinges, [output_dir filesep 'hinge_plots'] );
        [ hinge ] = fn_accept_hinge( element, hinge );
    end
end

%% Save Data
save([output_dir filesep 'model_analysis.mat'],'model')
save([output_dir filesep 'story_analysis.mat'],'story')
save([output_dir filesep 'element_analysis.mat'],'element')
save([output_dir filesep 'element_TH.mat'],'element_TH')
save([output_dir filesep 'element_PM.mat'],'element_PM')
save([output_dir filesep 'hinge_analysis.mat'],'hinge')
% writetable(element,[output_dir filesep 'element_linear.csv'])
end

