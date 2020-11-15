%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
analysis.model_id = 9;
analysis.proceedure = 'NDP';
analysis.id = 'IDA';
analysis.gm_set = 'FEMA_far_field';
analysis.run_z_motion = 1;

% Analysis options
analysis.summit = 0;
analysis.run_parallel = 0;
analysis.run_ida = 1;
analysis.post_process_ida = 1;
analysis.create_fragilities = 0;
analysis.plot_ida = 0;
analysis.detialed_post_process = 0;
analysis.run_sa_stripes = 1;
analysis.scale_increment = 0.25;
analysis.sa_stripes = [0.1 0.2 0.3 0.4 0.5];
analysis.collapse_drift = 0.10;
analysis.clear_existing_data = 0;
analysis.general_ida = 0;

% Secondary options
analysis.dead_load = 1;
analysis.live_load = 1;
analysis.opensees_SP = 1;
analysis.type = 1;
analysis.nonlinear = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.03;
analysis.hinge_stiff_mod = 10;
analysis.run_eigen = 0;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 1;
analysis.play_movie = 0;
analysis.movie_scale = 0;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';
analysis.joint_explicit = 0;
analysis.stories_nonlinear = inf;
analysis.stories_nonlinear_low = 0;
analysis.elastic_beams = 0;
analysis.simple_recorders = 0;

%% Import Packages
import ida.*

%% P695 Factors
% Max Dir Spectra (from Russ)
max_dir_spectra = readtable('+ida/ida_max_dir_spectra.csv');

% Period Based Ductility
mu_t_ew = 0.017 / 0.005; % rough estimate from pushover (may be slightly higher)
mu_t_ns = 0.004 / 0.001; % rough estimate from pushover (may be slightly higher)

% MCE
ida_results.direction = {'EW'; 'NS'};
ida_results.period = [1.14; 0.35];
ida_results.spectra = [0.56; 1.18]; % max direction spectra of ICSB motion
ida_results.mce = [0.79; 1.55]; % MCE Max from SP3, fixed to this site and model period

% SSF (based on table 7-1b)
SSF_ew = interp1([1.0 1.1 1.5 2 3 4 6 8], [1.00 1.07 1.14 1.19 1.27 1.32 1.41 1.49], mu_t_ew);
SSF_ns = interp1([1.0 1.1 1.5 2 3 4 6 8], [1.00 1.05 1.10 1.13 1.18 1.22 1.28 1.33], mu_t_ns);

%% Initial Setup
% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Define read and write directories
main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
model_dir = [main_dir '/' 'model_data'];
tcl_dir = [main_dir '/' 'opensees_data'];
asce41_dir = [main_dir '/' 'asce_41_data'];

% Load in Model Tables
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
load([asce41_dir filesep 'element_analysis.mat']);
load([asce41_dir filesep 'joint_analysis.mat']);

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);

%% Run Opensees Models
if analysis.run_ida || analysis.post_process_ida
    fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
end

%% Create Response and Consequence Fragilities
if analysis.create_fragilities
    write_dir = [main_dir '/' 'IDA' '/' 'Fragility Data'];
    if ~exist(write_dir,'dir')
        mkdir(write_dir)
    end
    fn_collect_ida_data(analysis, model, gm_set_table, ida_results, main_dir, write_dir)
    fn_create_fragilities(analysis, write_dir)
%     fn_LR_classification(analysis, model, gm_set_table)
end

%% Plot Results
if analysis.plot_ida
    fn_plot_ida(analysis, model, gm_set_table, ida_results, SSF_ew, main_dir)
    fn_plot_new_frags(analysis, model, main_dir)
%     fn_plot_hinge_values(analysis, model, IDA_scale_factors, gm_set_table)
end

%% Post Process Details of Specific GMs
if analysis.detialed_post_process
    % Inputs
    analysis.single_gm_scale = 2.9717;
    analysis.single_gm_set_id = 13;
    analysis.single_gm_pair_id = 2;
    
    % Get data
    fn_postprocess_single_IDA(analysis, model, gm_set_table)
end


