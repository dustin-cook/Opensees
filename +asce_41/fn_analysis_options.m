function [ analysis ] = fn_analysis_options( analysis )
% Description: Defaults secondary analysis options for ASCE 41 Assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs: Analysis Data Structure

% Outputs: Analysis Data Structure

% Assumptions:

%% Basic Defaults
% Run Options
analysis.run_opensees = 1; % 1 = Run opensees, 0 = use existing results
analysis.asce_41_post_process = 1; % 1 = run asce 41 post process logic
analysis.summit_SP = 0; % Write tcl files to be run on summit using OpenseesSP
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% Model Options
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)

% Opensees Analysis Options
analysis.ground_motion_scale_factor = 1; % Scale the GM amplitude
analysis.damping = 'rayleigh'; % rayleigh, modal, or simple
analysis.damp_ratio = 0.05; % Analysis damping ration
analysis.hinge_stiff_mod = 10; % Scale up stiffnes of hinges for a lumped plasticiy model. n value from Ibarra paper.
analysis.run_eigen = 1; % Run the eignen anlayis to get mode shapes and periods for the opensees analysis
analysis.initial_timestep_factor = 1; % reduction from eq timestep to analysis timestep
analysis.solution_algorithm = 1; % Run the opensees solution algorthm which will try different things 
analysis.collapse_drift = 0.20; % stop analysis at this drift and say collapse
analysis.joint_model = 2; % 1 = elastic elements, 2 = joint 3D
analysis.write_xml = 1; % Write and read opensees out files as xml files (0 = .txt files)
analysis.pushover_drift = 0.05; % Drfit limit where the pushover will go till
analysis.pushover_num_steps = 100; % Number of steps a pushover will take to get to the dirft limit

% Visuals and Graphics
analysis.element_plots = 0; % Plot hinge backnones and other per element visualizations
analysis.plot_recordings = 0; % Plot analysis results v recorded results
analysis.play_movie = 1; % Have opensees display a real time graphic of the building and analysis
analysis.movie_scale = 1; % Visual scale of the movie playback

%% Define Proceedure Options
if strcmp(analysis.proceedure,'test')
    analysis.type_list = 1; % just one linear dynamic analysis
    analysis.nonlinear_list = 0;
    analysis.dead_load_list = 1;
    analysis.live_load_list = 1;
    analysis.case_list = {'NA'};
elseif strcmp(analysis.proceedure,'NDP')
    analysis.type_list = [2, 2, 2, 1]; % Linear Pushover then NL Pushover x 2 then 1 NL dynamic
    analysis.nonlinear_list = [0, 1, 1, 1];
    analysis.dead_load_list = [1, 1, 1, 1];
    analysis.live_load_list = [1, 1, 1, 1];
    analysis.case_list = {'NA', 'NA', 'NA', 'NA'};
elseif strcmp(analysis.proceedure,'LDP') % Linear Test
    analysis.type_list = [1, 1]; % 1 = dynamic, 2 = pushover % 3 = static cyclic
    analysis.nonlinear_list = [0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.dead_load_list = [1.1, 0.9]; % Dead load factor for linear analysis
    analysis.live_load_list = [1.1, 0.0]; % Live load factor for linear analysis
    analysis.case_list = {'load_case_1', 'load_case_2'};
end


end

