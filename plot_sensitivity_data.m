clear all
close all
clc

% Collect results from model sensitivity study

% Define Model
analysis.model_id = 18;
analysis.proceedure = 'NDP';
analysis.id = 'baseline';
analysis.model = 'ICBS_model_5ew_col_base';

%% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Define outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study'];
plot_dir = [outputs_dir filesep 'summary_plots'];

%% Define sensitivity study parameters
model_name = {'both'};
x_var = {'energy_cov'};
x_lab = {'Strength & Deformation Capacity COV'};
% model_name = {'ductility', 'ductility_no_var', 'strength', 'both', 'both_more', 'num_bays'};
% x_var = {'ductility_cov', 'variant', 'strength_cov', 'energy_cov', 'ductility_cov', 'num_bays'};
% x_lab = {'Variation in Deformation Capacity', 'Level of Ductility', 'Variation in Strength','Variation in Strength & Deformation Capacity','Variation in Strength & Deformation Capacity', 'Number of Bays'};
% model_name = {'ductility', 'strength', 'both', 'num_bays'};
% x_var = {'ductility_cov', 'strength_cov', 'energy_cov', 'num_components'};
% x_lab = {'Deformation Capacity COV', 'Strength COV','Strength & Deformation Capacity COV', 'Number of Components'};
% y_var = {'med_sa_cp', 'med_sa_cp_10', 'med_sa_cp_25', 'med_sa_cp_50', 'med_sa_cp_75','med_sa_b', 'med_sa_drift_2', 'med_sa_drift_3', 'med_sa_drift_4', 'med_sa_lat_50', 'med_sa_lat_60', 'med_sa_lat_70', 'med_sa_lat_80', 'med_sa_lat_90', 'med_sa_grav_30', 'med_sa_grav_40', 'med_sa_grav_50', 'med_sa_grav_60', 'med_sa_grav_70', 'med_sa_grav_80'};

% model_name = { 'num_bays'};
% x_var = {'num_components'};
% x_lab = { 'Number of Components'};

% y_var = {'med_sa_grav_80','med_sa_drift_3','med_sa_lat_50'};
% y_var = {'med_sa_drift_mean'};
% y_var = {'med_sa_drift_mean_grav'};
% y_var = {'med_sa_drift_mean_grav_min'};
y_var = {'med_sa_cp'};

% y_var = {'med_sa_cp', 'med_sa_cp_10', 'med_sa_cp_25', 'med_sa_cp_50', 'med_sa_cp_75','med_sa_b', 'med_sa_drift_2', 'med_sa_drift_3', 'med_sa_drift_4', 'med_sa_lat_50', 'med_sa_lat_60', 'med_sa_lat_70', 'med_sa_lat_80', 'med_sa_lat_90', 'med_sa_grav_30', 'med_sa_grav_40', 'med_sa_grav_50', 'med_sa_grav_60', 'med_sa_grav_70', 'med_sa_grav_80','med_sa_drift_mean', 'med_sa_drift_mean_grav', 'med_sa_drift_mean_grav_min'};

%% Load sensitivity summary results
summary_results = readtable([outputs_dir filesep 'summary_results.csv']);

%% Collect data for each sensitivity study
for m = 1:length(model_name)
    % merge baselines data and model data
    base_data = summary_results(strcmp(summary_results.group,'baseline'),:);
    ICSB_2D = summary_results(strcmp(summary_results.group,'ICSB_2D'),:);
    ICSB_3D = summary_results(strcmp(summary_results.group,'ICSB_3D'),:);
    model_data = summary_results(strcmp(summary_results.group,model_name{m}),:);
    data = [base_data;  model_data];
%     data = summary_results;
    
    % Fit Linear Model
    for y = 1:length(y_var)
        data.ydata = data.med_sa_collapse_grav ./ data.(y_var{y});
        base_data.ydata = base_data.med_sa_collapse_grav ./ base_data.(y_var{y});
        ICSB_2D.ydata_grav = ICSB_2D.med_sa_collapse_grav ./ ICSB_2D.(y_var{y});
        ICSB_3D.ydata_grav = ICSB_3D.med_sa_collapse_grav ./ ICSB_3D.(y_var{y});
        ICSB_2D.ydata = ICSB_2D.med_sa_collapse ./ ICSB_2D.(y_var{y});
        ICSB_3D.ydata = ICSB_3D.med_sa_collapse ./ ICSB_3D.(y_var{y});
        mdl = fitlm(data, ['ydata~' x_var{m}]);

        % Plot trend with collapse margin
        mdl_plt = plot(mdl);
    %     scatter(data.(x_var{m}), data.collapse_margin_grav,'filled')
        plot_name = [model_name{m} ' - ' y_var{y}];
        xlabel(x_lab{m})
        ylabel('Collapse Indicator Ratio')
        title('') % delete title 
        delete(mdl_plt(3)) % delete CI lower
        delete(mdl_plt(4)) % delete CI upper
        legend('Sensitivity Models','Linear Trend');
        ylim([0.8, 2.2])
        hold on 
        scatter(base_data.(x_var{m}), base_data.ydata,100,'s','displayname','Baseline')
%         scatter(ICSB_2D.(x_var{m}), ICSB_2D.ydata,60,'o','m','filled','displayname','ICSB 2D Sidesway')
%         scatter(ICSB_2D.(x_var{m}), ICSB_2D.ydata_grav,60,'o','k','filled','displayname','ICSB 2D Gravity')
%         scatter(ICSB_3D.(x_var{m}), ICSB_3D.ydata,60,'d','m','filled','displayname','ICSB 3D Sideway')
%         scatter(ICSB_3D.(x_var{m}), ICSB_3D.ydata_grav,60,'d','k','filled','displayname','ICSB 3D Gravity')
%         legend off
        fn_format_and_save_plot( plot_dir, plot_name, 2 )
    end
end

