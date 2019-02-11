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
import asce_41.main_combine_load_case

%% Begin Method
% Pull in database of available models
if analysis.model_type == 1 % SDOF
    model_table = readtable(['inputs' filesep 'sdof_models.csv'],'ReadVariableNames',true);
elseif analysis.model_type == 2 % MDOF
    model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
end

% Pull in Element Database
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% Select Model for the analysis 
model = model_table(model_table.id == analysis.model_id,:);

% Create Analysis Directory
analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure];
if analysis.run_opensees && ~analysis.skip_2_outputs % Don't clear the file if you don't want to run opensees
    fn_make_directory( analysis.out_dir )
end

if ~analysis.skip_2_outputs % Don't skip to plotters
    % Run through all the steps of the procedure
    for i = 1:length(analysis.type_list)
        analysis.type = analysis.type_list(i);
        analysis.nonlinear = analysis.nonlinear_list(i);
        analysis.dead_load = analysis.dead_load_list(i);
        analysis.live_load = analysis.live_load_list(i);
        analysis.case = analysis.case_list{i};
        analysis.pushover_drift_x = analysis.pushover_drift_list_x(i);
        analysis.pushover_drift_z = analysis.pushover_drift_list_z(i);
        disp(['Running ' analysis.proceedure ' step ' num2str(i) ' of ' num2str(length(analysis.type_list)) ' ...'])

        %% Build Model
        disp('Building Model ...')
        main_build_model( model, analysis, ele_prop_table )

        %% Run and Postprocess Opensees Analysis
        disp('Running Opensees ...')
        main_opensees_analysis( model, analysis )

        %% Postprocess ASCE 41 data
        disp('Post Processing Via ASCE 41 ...')
        [ capacity(:,i) ] = main_ASCE_41_post_process( analysis, ele_prop_table );

        %% Analysis Checks
        disp('Validating Analysis Results ...')
        main_check_analysis( analysis, ele_prop_table, capacity, i )
    end

    %% Combine Load Cases
    if strcmp(analysis.proceedure,'LDP')
        main_combine_load_case( analysis )
    end
end

%% Compile Results and Create Visuals
disp('Plotting Analysis Results ...')
main_plot_analysis_results( model, analysis, ele_prop_table )

%% LaTeX Report Writer


disp('Analysis Complete!')
end

