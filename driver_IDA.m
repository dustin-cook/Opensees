%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% Initial Setup
% Define Model
analysis.model_id = 6;
analysis.proceedure = 'NDP';
analysis.summit = 1;
analysis.run_ida = 1;
analysis.gm_set = 'FEMA_far_field';

% IDA Inputs
hazard.curve.rp = [43, 72, 224, 475, 975, 2475, 4975];
hazard.curve.pga = [0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];
analysis.collapse_drift = 0.1;

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

% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Import packages
import plotting_tools.fn_format_and_save_plot

% Set up Parallel Workers
parpool;

%% Define read and write directories
model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'model_data'];
tcl_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'opensees_data'];

% Load in Model Tables
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
element = readtable([model_dir filesep 'element.csv'],'readVariableNames',true);

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);
IDA_scale_factors = hazard.curve.pga ./ gm_median_pga;

%% Run Opensees Model
if analysis.run_ida
    for i = 1:length(IDA_scale_factors)
        scale_factor = IDA_scale_factors(i);
        for gms = 1:height(gm_set_table)
            fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, tcl_dir)
            clc 
        end
    end
end

%% Plot Results
plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'IDA' '/' 'IDA Plots'];
id = 0;
for i = 1:length(IDA_scale_factors)
    for gms = 1:height(gm_set_table)
        for d = 1:2
            id = id + 1;
            % Load data
            outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '/' 'IDA' '/' 'Scale_' num2str(IDA_scale_factors(i)) '/' 'GM_' num2str(gm_set_table.set_id(gms)) '_' num2str(d)];
            load([outputs_dir filesep 'summary_results.mat'])
            ida.id(id,1) = id;
            ida.scale(id,1) = IDA_scale_factors(i);
            ida.gm_name{id,1} = gm_set_table.eq_name{gms};
            ida.sa_x(id,1) = summary.sa_x;
            ida.sa_z(id,1) = summary.sa_z;
            ida.drift_x(id,1) = summary.max_drift_x;
            ida.drift_z(id,1) = summary.max_drift_z;
            ida.collapse(id,1) = summary.collapse;
        end
    end
end

% End Parallel Process
delete(gcp('nocreate'))

% Plot IDA curves
hold on
for gms = 1:height(gm_set_table)
    plot(ida.drift_x(strcmp(ida.gm_name,gm_set_table.eq_name{gms})),ida.sa_x(strcmp(ida.gm_name,gm_set_table.eq_name{gms})))
end
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( plot_dir, 'IDA Plot EW Frame Direction', 2 )
hold on
for gms = 1:height(gm_set_table)
    plot(ida.drift_z(strcmp(ida.gm_name,gm_set_table.eq_name{gms})),ida.sa_z(strcmp(ida.gm_name,gm_set_table.eq_name{gms})))
end
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( plot_dir, 'IDA Plot NS Wall Direction', 2 )

% Save results as csv
ida_table = struct2table(ida);
writetable(ida_table,[plot_dir filesep 'ida_results.csv'])

