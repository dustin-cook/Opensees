%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
% analysis.model_id = 4;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'P58'; 
analysis.id = 'Dissertation_Study_LP'; % ID of the analysis for it to create its own directory
% analysis.gm_set = 'FEMA_far_field';
analysis.run_z_motion = 0;

% Analysis options
% analysis.summit = 0;
% analysis.run_parallel = 0;
% analysis.run_ida = 0;
% analysis.post_process_ida = 1;
% analysis.create_fragilities = 0;
% analysis.plot_ida = 0;
% analysis.detialed_post_process = 0;
% analysis.run_sa_stripes = 1;
% analysis.scale_method = '2D'; % 'maxdir' or 'geomean' oe 2D
% analysis.scale_increment = 0.25;
% analysis.sa_stripes = [0.2 0.4];
% analysis.collapse_drift = 0.10; % Drift ratio at which we are calling the model essentially collapsed
% analysis.clear_existing_data = 0;
% analysis.general_ida = 0;
% analysis.write_xml = 1;

% Nonlinear options
analysis.nonlinear = 1;
analysis.nonlinear_type = 'lumped'; % lumped or fiber
analysis.stories_nonlinear = inf;
analysis.stories_nonlinear_low = 0;
analysis.elastic_beams = 0;
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;

% Secondary options
analysis.dead_load = 1; % Dead load factor
analysis.live_load = 0.2; % live load factor
% analysis.opensees_SP = 0;
analysis.type = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.025;
% analysis.solution_algorithm = 1;
% analysis.initial_timestep_factor = 1;
% analysis.suppress_outputs = 1;
% analysis.play_movie = 0;
% analysis.movie_scale = 10;
% analysis.algorithm = 'KrylovNewton';
% analysis.integrator = 'Newmark 0.5 0.25';
analysis.joint_model = 1;
analysis.joint_explicit = 0;
analysis.simple_recorders = 1;
analysis.additional_elements = 1;

analysis.build_model_remote = 1;


% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(1:24,:);
% model_data = model_data(1:4,:);
% model_data = model_data(model_data.num_stories == 4,:);
% model_data = model_data(model_data.ie == 1,:);
% model_data = model_data(model_data.design_drift == 0.02,:);
num_models = height(model_data);

% Define remote directory
remote_dir = ['C:\Users\dtc2\OneDrive - NIST\NIST\Dissertation Archetype Study\Models\Archetypes RCMF\EDPs_2D'];

%% Import Packages
import build_model.main_build_model

%% Go through each model and build model files
for m = 1:num_models % run for each model     
    % Set basic model data
    model = model_data(m,:);
    analysis.model_id = model.id;
    fprintf('Building Model %i of %i: %s\n', m, num_models, model.name{1})
    
    % Define Model Directory
    if analysis.build_model_remote
        % Define Remote Directory
        analysis.out_dir = [remote_dir filesep model.id{1}];
%         if ~exist(model_remote_dir,'dir')
%             mkdir(model_remote_dir)
%         end
    else
        analysis.out_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
    end

    % Create model tables
    main_build_model( model, analysis, [] )
end
