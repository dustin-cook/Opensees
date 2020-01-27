clear all
close all
clc

% Collect results from model sensitivity study

% Define Model
analysis.model_id = 18;
analysis.proceedure = 'NDP';
analysis.id = 'baseline_beams_include';
analysis.model = 'ICBS_model_5ew_col_base';

%% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Define outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study'];
plot_dir = [outputs_dir filesep 'summary_plots'];

%% Define sensitivity study parameters
model_name = {'both'};
x_var = {'ductility_cov'};
x_lab = {'Variation in Strength & Ductility'};
% model_name = {'ductility', 'ductility_no_var', 'strength', 'both', 'both_more', 'num_bays'};
% x_var = {'ductility_cov', 'variant', 'strength_cov', 'ductility_cov', 'ductility_cov', 'num_bays'};
% x_lab = {'Variation in Ductility', 'Level of Ductility', 'Variation in Strength','Variation in Strength & Ductility','Variation in Strength & Ductility', 'Number of Bays'};
% model_name = {'ductility', 'strength', 'both', 'num_bays'};
% x_var = {'ductility_cov', 'strength_cov', 'ductility_cov', 'num_bays'};
% x_lab = {'Variation in Ductility', 'Variation in Strength','Variation in Strength & Ductility', 'Number of Bays'};

%% Load sensitivity summary results
summary_results = readtable([outputs_dir filesep 'summary_results.csv']);

%% Collect data for each sensitivity study
for m = 1:length(model_name)
    % merge baselines data and model data
    base_data = summary_results(strcmp(summary_results.group,'baseline'),:);
    model_data = summary_results(strcmp(summary_results.group,model_name{m}),:);
    data = [base_data;  model_data];
    
    % Fit Linear Model
    mdl = fitlm(data, ['collapse_margin_grav~' x_var{m}]);
    
    % Plot trend with collapse margin
    plot(mdl)
%     scatter(data.(x_var{m}), data.collapse_margin_grav,'filled')
    plot_name = model_name{m};
    xlabel(x_lab{m})
    ylabel('Ratio of Collapse to CP')
    title('')
    legend('Data',['Fit: R-squared = ', num2str(round(mdl.Rsquared.ordinary,2))],'95% Confidence Bounds');
    fn_format_and_save_plot( plot_dir, plot_name, 2 )
end

