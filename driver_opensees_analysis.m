% Run Tcl file in opensees
clear
close
clc

%% DEFINE INPTUTS
% Primary Inputs
analysis.model_id = 1;
analysis.gm_seq_id = 6;
analysis.name = 'test_analysis_1';

% Secondary Inputs
analysis.type = 3;
analysis.model_type = 1;
analysis.max_disp = 7;
analysis.time_step = 0.01;
analysis.nonlinear = 0;
analysis.dead_load = 1.0;
analysis.live_load = 1.0;
analysis.accidental_torsion = 0;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.05;
analysis.hinge_stiff_mod = 10;
analysis.play_movie = 1;
analysis.run_eigen = 1;
analysis.run_opensees = 1;
analysis.initial_timestep_factor = 1;
analysis.solution_algorithm = 0;

%% Initial Setup
import opensees.main_opensees_analysis

%% Run Opensees Analysis
main_opensees_analysis( analysis )




