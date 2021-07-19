function [ ] = main_ASCE_7( analysis, load_case_id, model )
% Description: Main script facilitating an ASCE 41 teir 3 assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
import opensees.main_opensees_analysis
import asce_7.fn_combine_load_cases

% Create Analysis Directory
analysis.model_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep 'model_data'];
analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep load_case_id];
fn_make_directory( analysis.out_dir )

%% Run each load combo
for i = 1:length(analysis.type_list)
    analysis.type = analysis.type_list(i);
    analysis.nonlinear = analysis.nonlinear_list(i);
    analysis.run_drifts = analysis.drift_run_list(i);
    analysis.dead_load = analysis.dead_load_list(i);
    analysis.live_out_load = analysis.live_out_load_list(i);
    analysis.live_in_load = analysis.live_in_load_list(i);
    analysis.eq_lat_load_factor = analysis.eq_lat_load_list(i);
    analysis.eq_vert_load_factor = analysis.eq_vert_load_list(i);
    disp(['Running ' analysis.proceedure ' step ' num2str(i) ' of ' num2str(length(analysis.type_list)) ' ...'])

    % Run and Postprocess Opensees Analysis
    disp('Running Opensees ...')
    main_opensees_analysis( model, analysis )

    % Combine Load cases
    disp('Combing load cases ...')
    fn_combine_load_cases( analysis, i );
end

disp('Analysis Complete!')
end

