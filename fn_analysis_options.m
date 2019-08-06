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
analysis.run_opensees_post_process = 1; % 1 = Run opensees, 0 = use existing results
analysis.asce_41_post_process = 1; % 1 = run asce 41 post process logic
analysis.opensees_SP = 1; % 0 = Standard OpenSees; 1 = OpenseesSP

% Model Options
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0;

% Opensees Analysis Options
analysis.ground_motion_scale_factor = 1.167259; % Scale the GM amplitude   2.9717
analysis.damping = 'rayleigh'; % rayleigh, modal, or simple
analysis.hinge_stiff_mod = 10; % Scale up stiffnes of hinges for a lumped plasticiy model. n value from Ibarra paper.
analysis.run_eigen = 1; % Run the eignen anlayis to get mode shapes and periods for the opensees analysis
analysis.initial_timestep_factor = 1; % reduction from eq timestep to analysis timestep
analysis.solution_algorithm = 1; % Run the opensees solution algorthm which will try different things 
analysis.collapse_drift = 0.06; % stop analysis at this drift and say collapse
analysis.joint_model = 1; % 1 = beam/column elements, 2 = joint 3D
analysis.joint_explicit = 0; % 0 = rigid, 1 = model joint nonlinearity (could automate this based on first assessment of joints)
analysis.write_xml = 1; % Write and read opensees out files as xml files (0 = .txt files, which is currently broken, on purpose)
analysis.pushover_num_steps = 200; % Number of steps a pushover will take to get to the dirft limit
analysis.cyclic_pushover_peak_drifts = [0.4, 0.5, 0.6]; % Percent of the final Pushover drift of each cycle
analysis.hinge_group_length = 10;
analysis.filter_accel = 0;
analysis.filter_freq_range = [0.5, 1.5];
analysis.simple_recorders = 0;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

% Visuals and Graphics
analysis.element_plots = 1; % Plot hinge backnones and other per element visualizations
analysis.plot_recordings = 1; % Plot analysis results v recorded results
analysis.play_movie = 1; % Have opensees display a real time graphic of the building and analysis
analysis.movie_scale = 1; % Visual scale of the movie playback
analysis.hinge_stories_2_plot = 1;
analysis.suppress_outputs = 0;

%% Define Proceedure Options
if strcmp(analysis.proceedure,'test')
    analysis.type_list = [1, 1];
    analysis.nonlinear_list = [1, 1];
    analysis.dead_load_list = [1, 1];
    analysis.live_load_list = [1, 1];
    analysis.case_list = {'backbones_pushover', 'backbones'};
    analysis.pushover_drift_list_x = [0.01, 0.01];
    analysis.pushover_drift_list_z = [0.05, 0.05];
    analysis.accidental_torsion_list = [0, 0];
    analysis.damp_ratio_list = [0.03, 0.03]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'torsion')
    analysis.type_list = [1, 1];
    analysis.nonlinear_list = [0, 0];
    analysis.dead_load_list = [1, 1];
    analysis.live_load_list = [1, 1];
    analysis.case_list = {'NA', 'torsion_check'};
    analysis.pushover_drift_list_x = [NaN, NaN];
    analysis.pushover_drift_list_z = [NaN, NaN];
    analysis.accidental_torsion_list = [0, 1];
    analysis.damp_ratio_list = [0.028, 0.028];
    
elseif strcmp(analysis.proceedure,'Pushover')
    analysis.type_list = [2, 2]; % NonLinear Pushover then NL Pushover x 2
    analysis.nonlinear_list = [1, 1];
    analysis.dead_load_list = [1, 1];
    analysis.live_load_list = [1, 1];
    analysis.case_list = {'backbones_pushover', 'backbones'};
    analysis.pushover_drift_list_x = [0.05, 0.05]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.05 0.05];
    analysis.accidental_torsion_list = [0, 0];
    analysis.damp_ratio_list = [0.03, 0.03]; % Analysis damping ratio
    
elseif strcmp(analysis.proceedure,'NDP')
    analysis.type_list = [2, 2, 1]; % Linear Pushover then NL Pushover x 2 then 1 NL dynamic
    analysis.nonlinear_list = [1, 1, 1];
    analysis.dead_load_list = [1, 1, 1];
    analysis.live_load_list = [1, 1, 1];
    analysis.case_list = {'backbones_pushover', 'backbones', 'NA'};
    analysis.pushover_drift_list_x = [0.02, 0.02, NaN]; % Drift limit where the pushover will go till
    analysis.pushover_drift_list_z = [0.01, 0.01, NaN];
    analysis.accidental_torsion_list = [0, 0, 1];
    analysis.damp_ratio_list = [0.03, 0.03, 0.03]; % Analysis damping ratio
   
elseif strcmp(analysis.proceedure,'LDP') % Linear Test
    analysis.type_list = [1, 1]; % 1 = dynamic, 2 = pushover % 3 = static cyclic
    analysis.nonlinear_list = [0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
    analysis.dead_load_list = [1.1, 0.9]; % Dead load factor for linear analysis
    analysis.live_load_list = [1.1, 0.0]; % Live load factor for linear analysis
    analysis.case_list = {'load_case_1', 'load_case_2'};
    analysis.pushover_drift_list = [NaN, NaN];
    analysis.accidental_torsion_list = [0, 0];
    analysis.damp_ratio = [0.05, 0.05]; % Analysis damping ratio
end


end

