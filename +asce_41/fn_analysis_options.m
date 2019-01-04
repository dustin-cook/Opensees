function [ analysis ] = fn_analysis_options( analysis )
% Description: Defaults secondary analysis options for ASCE 41 Assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs: Analysis Data Structure

% Outputs: Analysis Data Structure

% Assumptions:

%% Basic Defaults
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)
analysis.ground_motion_scale_factor = 1; % Scale the GM amplitude
analysis.damping = 'rayleigh'; % rayleigh, modal, or simple
analysis.damp_ratio = 0.05;
analysis.hinge_stiff_mod = 10; % Scale up stiffnes of hinges for a lumped plasticiy model. n value from Ibarra paper.
analysis.play_movie = 1;
analysis.movie_scale = 1; % Visual scale of the movie playback
analysis.run_eigen = 1;
analysis.initial_timestep_factor = 1; % reduction from eq timestep to analysis timestep
analysis.solution_algorithm = 1;
analysis.collapse_drift = 0.20; % stop analysis at this drift and say collapse
analysis.joint_model = 2; % 1 = elastic elements, 2 = joint 3D
analysis.full_recorders = 0; % 0 = simple recorders, 1 = full recorders
analysis.write_xml = 1; % Write and read opensees out files as xml files (0 = .txt files)
analysis.plot_hinges = 1; % Plot hinge backnones
analysis.plot_asce = 1;


%% Define Proceedure Options
if strcmp(analysis.proceedure,'LDP')
    analysis.type_list = 1; % 1 = dynamic, 2 = pushover % 3 = static cyclic
    analysis.nonlinear_list = 0; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
elseif strcmp(analysis.proceedure,'NDP')
    analysis.type_list = [2, 2, 2, 1]; % Pushover x 2 then 1 dynamic
    analysis.nonlinear_list = [0, 1, 1, 1];
elseif strcmp(analysis.proceedure,'test') % Linear Test
    analysis.type_list = 1; % Pushover x 2 then 1 dynamic
    analysis.nonlinear_list = 0;
end

end

