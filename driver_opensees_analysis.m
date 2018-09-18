% Run Tcl file in opensees
clear
close
clc

%% DEFINE INPTUTS
% Primary Inputs
analysis.model_id = 10;
analysis.gm_seq_id = 8;
analysis.name = 'test';

% Secondary Inputs
analysis.type = 2; % 1 = dynamic, 2 = pushover
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF
analysis.pushover_drift = 0.10;
analysis.pushover_num_steps = 100;
analysis.pushover_direction = 'x';
analysis.time_step = 0.01;
analysis.nonlinear = 0;
analysis.dead_load = 1.0;
analysis.live_load = 1.0;
analysis.accidental_torsion = 0;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.05;
analysis.hinge_stiff_mod = 10;
analysis.play_movie = 1;
analysis.movie_scale = 1;
analysis.run_eigen = 1;
analysis.run_opensees = 1;
analysis.initial_timestep_factor = 1;
analysis.solution_algorithm = 0;
analysis.joint_model = 2; % 1 = elastic elements, 2 = joint 3D
analysis.full_recorders = 0; % 0 = simple recorders, 1 = full recorders

%% Initial Setup
import opensees.main_opensees_analysis

%% Run Opensees Analysis
main_opensees_analysis( analysis )




