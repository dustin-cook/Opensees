function [ ] = main_ASCE_41( analysis )
% Description: Main script facilitating an ASCE 41 teir 3 assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
import build_model.main_build_model
import opensees.main_opensees_analysis
import asce_41.main_ASCE_41_post_process
import asce_41.main_plot_analysis_results


%% Begin Method
% Pull in database of available models
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);

% Select Model for the analysis 
model = model_table(model_table.id == analysis.model_id,:);

% Create Outputs Directory
output_dir = ['outputs/' model.name{1} filesep 'model data'];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Run through all the steps of the procedure
for i = 1:length(analysis.type_list)
    analysis.type = analysis.type_list(i);
    analysis.nonlinear = analysis.nonlinear_list(i);
    
    %% Build Model
    main_build_model( model, analysis, output_dir )
    
    %% Run and Postprocess Opensees Analysis
    main_opensees_analysis( analysis )

    %% Postprocess ASCE 41 data
    main_ASCE_41_post_process( analysis )

    %% Analysis Checks
    main_check_analysis( )
end

%% Compile Results and Create Visuals
main_plot_analysis_results( analysis )
    
%% Latex Report Writer


end

