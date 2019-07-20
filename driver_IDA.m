%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
analysis.model_id = 11;
analysis.proceedure = 'NDP';
analysis.id = 2;
analysis.summit = 0;
analysis.run_ida = 0;
analysis.post_process_ida = 1;
analysis.plot_ida = 0;
analysis.gm_set = 'FEMA_far_field';

% IDA Inputs
% hazard.curve.rp = [22, 35];%, 64, 108, 144];
% hazard.curve.pga = [0.128, 0.192];%, 0.288, 0.376, 0.425];
hazard.curve.rp = [22, 35, 43, 64, 72, 108, 144, 224, 475, 975, 2475, 4975];
hazard.curve.pga = [0.128, 0.192, 0.224, 0.288, 0.308, 0.376, 0.425 0.502, 0.635, 0.766, 0.946, 1.082];
% hazard.curve.rp = [43, 72, 224, 475, 975, 2475, 4975];
% hazard.curve.pga = [0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];
analysis.collapse_drift = 0.06;

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
analysis.simple_recorders = 1;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 1;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

%% P695 Factors
% Max Dir Spectra (from Russ)
max_dir_spectra = readtable('ida_max_dir_spectra.csv');

% Period Based Ductility
mu_t_ew = 0.016 / 0.005; % rough estimate from pushover (may be slightly higher)
mu_t_ns = 0.004 / 0.0017; % rough estimate from pushover (may be slightly higher)

% MCE
ida_results.direction = {'EW'; 'NS'};
ida_results.period = [1.14; 0.44];
ida_results.spectra = [0.56; 0.79]; % max direction spectra
ida_results.mce = [0.53; 1.36]; % MCE Max from SP3, fixed to this site and model period

% SSF (based on table 7-1b)
SSF_ew = interp1([3,4], [1.275, 1.33],mu_t_ew);
SSF_ns = interp1([1.0 1.1 1.5 2 3 4 6 8], [1.00 1.05 1.1 1.13 1.18 1.22 1.28 1.33],mu_t_ew);

% Dispersion
beta_rtr = 0.4;
beta_tot = 0.6;

%% Initial Setup
% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Define read and write directories
model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'model_data'];
tcl_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'opensees_data'];
asce41_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'asce_41_data'];

% Load in Model Tables
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
load([asce41_dir filesep 'element_analysis.mat']);

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);
IDA_scale_factors = hazard.curve.pga ./ gm_median_pga;

%% Run Opensees Models
if analysis.run_ida || analysis.post_process_ida
% parpool; % Set up Parallel Workers
for i = 1:length(IDA_scale_factors)
    error_count = 0;
    scale_factor = IDA_scale_factors(i);
    analysis.ground_motion_scale_factor = scale_factor;
    run_ida = analysis.run_ida;
    for gms = 1:height(gm_set_table)
        % Run Opensees
        if run_ida
            % Suppress MATLAB warnings
            warning('off','all')
            fprintf('Running Scale Factor %4.2f for Ground Motion ID: %i-%i \n\n', scale_factor, gm_set_table.set_id(gms), gm_set_table.pair(gms))
            [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, ida_results.period, tcl_dir);
            if exit_status == 1
                error_count = error_count + 1;
            end
        else
            exit_status = 0;
        end

        if analysis.post_process_ida && exit_status ~= 1
%             try
                fprintf('Postprocessing Opensees Ouputs\n')
                fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor)
%             catch
%                 error_count = error_count + 1;
%             end
        end
        fprintf('\n')
    end
    fprintf('%i Failed GMs for Scale Factor %4.2f \n\n', error_count, scale_factor)
end
delete(gcp('nocreate')) % End Parallel Process
end

%% Plot Results
if analysis.plot_ida
    fn_plot_ida(analysis, model, IDA_scale_factors, gm_set_table, max_dir_spectra, ida_results, SSF_ew, SSF_ns)
end


