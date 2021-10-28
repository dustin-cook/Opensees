%% Clear the Workspace
clear
close
clc
fclose('all');

%% Description: Method to build an Opensees model and run a ASCE 41-17 teir 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% User Inputs (Think about changing this to a file read and command line execution)
% analysis.model_id = 2;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'Pushover'; % LDP or NDP or test
analysis.id = 'Dissertation_Study'; % ID of the analysis for it to create its own directory
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% dynamic analysis inputs
% analysis.gm_seq_id = 16; % Maybe also make this part ot the defaults or model?

% Archetpye Inputs
analysis.additional_elements = 1; % this is the leaning column
analysis.eq_lat_load_factor = 1;

%% Initial Setup
% Import packages
import asce_41.main_ASCE_41

% Define remote directory
remote_dir = ['G:\My Drive\Dissertation Archetype Study\Archetypes RCMF\Archetype Model Responses'];

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Pull Model Data
% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(model_data.num_stories == 12,:);
% model_data = model_data(model_data.num_stories ~= 20,:);
% model_data = model_data(model_data.ie ~= 1.25,:);
model_data = model_data(~contains(model_data.name,'drift'),:);
num_models = height(model_data);

%% Initiate Analysis
for m = 2:num_models % run for each model    
    % Load basic model data
    model = model_data(m,:);
    analysis.model_id = model.id;
    fprintf('Running Model %i of %i: %s\n', m, num_models, model.name{1})

    % Define read and write directories
%     main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
%     model_dir = [main_dir '/' 'model_data'];
%     if ~exist(model_dir,'dir')
%         mkdir(model_dir)
%     end
%     tcl_dir = [main_dir '/' 'opensees_data'];
%     if ~exist(tcl_dir,'dir')
%         mkdir(tcl_dir)
%     end
%     model_remote_dir = [remote_dir filesep model.name{1}];
%     if ~exist(model_remote_dir,'dir')
%         mkdir(model_remote_dir)
%     end
%     ELFP_model_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'model_data'];
    % tcl_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'eigen_analysis'];
    % asce41_dir = [main_dir '/' 'asce_41_data'];
    
    tic
    main_ASCE_41( analysis )
    toc
end

