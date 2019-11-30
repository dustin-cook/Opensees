%% Script to:
% 1. pull sensitivity results from Dropbox
% 2. run frag curve post processors on all of them
% 3. save results in repo folders
% 4. Create cummulative plots 

clear all
close all
clc

import ida.*
import plotting_tools.*

%% Define inputs
model_baseline_name = 'NDP_baseline_1';
model.name{1} = 'ICBS_model_5ew_col_base';
model_names = {'Model_uniform' 'NDP_baseline' 'Model' 'Model' 'Model' 'Model' 'Model' 'NDP_Model_extreme' };
model_ids = {'1' '1' '1' '4' '6' '8' '10' '1'};

% model_names = {'NDP_Model_extreme'};
% model_ids = {'1'};

analysis.run_z_motion = 0;

ida_results.direction = {'EW'; 'NS'};
ida_results.period = [1.14; 0.35];
ida_results.spectra = [0.56; 1.18]; % max direction spectra of ICSB motion
ida_results.mce = [0.79; 1.55]; % MCE Max from SP3, fixed to this site and model period

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep 'FEMA_far_field' filesep 'ground_motion_set.csv'],'ReadVariableNames',true);

%% For each model
for m = 1:length(model_names)
    % Define Analysis Name and Location
    analysis.proceedure = model_names{m};
    analysis.id = model_ids{m};
    write_dir = ['outputs/' model.name{1} '/' analysis.proceedure '_' analysis.id '/IDA/Fragility Data'];
    pushover_dir = ['outputs' '/' model.name{1} '/' model_baseline_name '/' 'pushover'];
    model_dir = ['outputs' '/' model.name{1} '/' model_baseline_name '/' 'opensees_data'];
    if ~exist(write_dir,'dir')
        mkdir(write_dir)
    end
    
    % Post Process for Fragilities
    sprintf('%s_%s', analysis.proceedure, analysis.id)
%     fn_collect_ida_data(analysis, model, gm_set_table, ida_results, write_dir, pushover_dir, model_dir)
%     fn_create_fragilities(analysis, gm_set_table, write_dir)
 
    % Collect Summary Statistics
    read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
    ele_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'asce_41_data'];
    [percent_CP, sa_med_col, sa_med_cp, col_margin, min_cp_med, mean_cp_med, max_cp_med, column_cov, column_base_cov, column_min, column_range] = fn_collect_cummary_data(read_dir, ele_dir);
    models.percent_CP(:,m) = percent_CP;
    models.sa_med_col(m) = sa_med_col;
    models.sa_med_cp(m) = sa_med_cp;
    models.col_margin(m) = col_margin;
    models.min_cp_med(m) = min_cp_med;
    models.mean_cp_med(m) = mean_cp_med;
    models.max_cp_med(m) = max_cp_med;
    models.column_cov(m) = column_cov;
    models.column_base_cov(m) = column_base_cov;
    models.column_min(m) = column_min;
    models.column_range(m) = column_range;
end

%% Create Plots
outputs_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep 'NDP_IDA_new' filesep 'sensitivity_study'];

rank = ((1:44)/44)';
cmap = colormap(jet(100));
hold on
for i = 1:length(model_names)
    plot(models.percent_CP(:,i),rank,'color',cmap(max(round(100*models.column_cov(i)),1),:),'linewidth',1.25)
end
xlabel('Fraction of Components Exceeding CP')
ylabel('P[Collapse]')
h = colorbar;
h.Limits = [0,0.7];
ylabel(h, 'Column Deformation Capacity COV')
grid on
box on
fn_format_and_save_plot( outputs_dir, 'Multi Percent CP plot', 4, 1 )

hold on
plot(models.column_cov, models.sa_med_cp,'b','marker','o','DisplayName','CP')
plot(models.column_cov, models.sa_med_col,'k','marker','d','DisplayName','Collapse')
ylim([0,1])
xlabel('Column Deformation Capacity COV')
ylabel('Median Sa')
legend('location','northeast')
grid on
box on
fn_format_and_save_plot( outputs_dir, 'Median Exceedence Trend', 4, 1 )

plot(models.column_cov, models.col_margin,'k','marker','o')
ylim([0,2])
xlabel('Column Deformation Capacity COV')
ylabel('Collapse Margin')
grid on
box on
fn_format_and_save_plot( outputs_dir, 'Collapse Margin Trend', 4, 1 )

hold on
plot(models.column_cov, models.min_cp_med,'b','marker','o','DisplayName','Min')
plot(models.column_cov, models.mean_cp_med,'k','marker','d','DisplayName','Mean')
plot(models.column_cov, models.max_cp_med,'r','marker','s','DisplayName','Max')
ylim([0,2])
xlabel('Column Deformation Capacity COV')
ylabel('Median Deformation Demand to CP Ratio')
legend('location','northeast')
grid on
box on
fn_format_and_save_plot( outputs_dir, 'Min Mean Max Trend', 4, 1 )

function [percent_CP, sa_med_col, sa_med_cp, col_margin, min_cp_med, mean_cp_med, max_cp_med, column_cov, column_base_cov, column_min, column_range] = fn_collect_cummary_data(read_dir, ele_dir)
load([read_dir filesep 'frag_curves.mat'])
load([read_dir filesep 'new_frag_curves.mat'])
load([read_dir filesep 'gm_data.mat'])
load([ele_dir filesep 'element_analysis.mat'])

% Percent Comps exceede CP discrete fragilitu
percent_CP = sort(gm_data.collapse.cols_walls_1_percent_cp);

% Median Collapse Capacity
sa_med_col = frag_curves.collapse.theta;

% Collapse Margin
sa_med_cp = frag_curves.cols_walls_1.cp.theta(1);
col_margin = sa_med_col / sa_med_cp;

% Min Mean Max CP
min_cp_med = new_frag_curves.collapse.cols_walls_1_min_cp.theta;
mean_cp_med = new_frag_curves.collapse.cols_walls_1_mean_cp.theta;
max_cp_med = new_frag_curves.collapse.cols_walls_1_max_cp.theta;

% Calculate Element COV
columns = element(element.story == 1 & strcmp(element.type,'column'),:);
column_cov = std(columns.b_hinge_1 + columns.b_hinge_2) / mean(columns.b_hinge_1 + columns.b_hinge_2);
column_base_cov = std(columns.b_hinge_1) / mean(columns.b_hinge_1);
column_min = min(columns.b_hinge_1);
column_range = max(columns.b_hinge_1) - min(columns.b_hinge_1);

end