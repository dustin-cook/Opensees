clear all
close all
clc

import asce_41.*
import plotting_tools.*
import opensees.write_tcl.*

% Define Model
analysis.model_id = 18;
analysis.proceedure = 'NDP';
analysis.id = 'baseline';
analysis.model = 'ICBS_model_5ew_col_base';

% Secondary options
analysis.type = 1;
analysis.nonlinear = 1;
analysis.stories_nonlinear = 1; % Default to all modeling all stories as nonlinear when doing NDP
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;
analysis.suppress_outputs = 1;
analysis.joint_explicit = 0;
analysis.joint_model = 1; % 1 = beam/column elements, 2 = joint 3D

%% Define inputs and outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
inputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'asce_41_data'];
os_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'opensees_data'];
os_model_csv_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'model_data'];
outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study'];
new_models_dir = [outputs_dir filesep 'Model Files'];
mkdir(outputs_dir)
mkdir(new_models_dir)

%% Define Element Properties
% Load Data
ele_prop_table = readtable(['inputs' filesep 'element.csv']);
load([inputs_dir filesep 'element_analysis.mat'])
load([inputs_dir filesep 'joint_analysis.mat'])
load([inputs_dir filesep 'hinge_analysis.mat'])
load([inputs_dir filesep 'node_analysis.mat'])
load([inputs_dir filesep 'story_analysis.mat'])

% Filter to just first story columns
first_story_columns = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,first_story_columns.node_1),:);

% %% Get column hinge properties
% for e = 1:height(first_story_columns)
%     for s = 1:2
%         % Define Model Parameters
%         ele = first_story_columns(e,:);
%         ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
%         first_story_columns.(['a_mean_' num2str(s)])(e,:) = ele.(['a_hinge_' num2str(s)]);
%         first_story_columns.(['b_mean_' num2str(s)])(e,:) = ele.(['b_hinge_' num2str(s)]);
%         first_story_columns.(['a_std_' num2str(s)])(e,:) = ele.(['a_hinge_' num2str(s)])*(1 - 0.5264);
%         first_story_columns.(['b_std_' num2str(s)])(e,:) = ele.(['b_hinge_' num2str(s)])*(1 - 0.5920);
%         first_story_columns.(['a_mu_' num2str(s)])(e,:)= log(ele.(['a_hinge_' num2str(s)]).^2./sqrt(first_story_columns.(['a_std_' num2str(s)])(e,:).^2 + ele.(['a_hinge_' num2str(s)]).^2));
%         first_story_columns.(['b_mu_' num2str(s)])(e,:) = log(ele.(['b_hinge_' num2str(s)]).^2./sqrt(first_story_columns.(['b_std_' num2str(s)])(e,:).^2 + ele.(['b_hinge_' num2str(s)]).^2));
%         first_story_columns.(['a_beta_' num2str(s)])(e,:) = sqrt(log(first_story_columns.(['a_std_' num2str(s)])(e,:).^2./ele.(['a_hinge_' num2str(s)]).^2+1));
%         first_story_columns.(['b_beta_' num2str(s)])(e,:) = sqrt(log(first_story_columns.(['b_std_' num2str(s)])(e,:).^2./ele.(['b_hinge_' num2str(s)]).^2+1));
%     end
% end
% 
% % Save Simulated Data
% save([outputs_dir filesep 'sim_data.mat'],'first_story_columns')
% 
% %% Simulate 1000 models
% % Simulate diversity of each model
params = {'a', 'b'};
% for sim = 1:1000
%     models.id(sim,1) = sim;
%     P = rand(height(first_story_columns),1);
%     for p = 1:length(params)
%         for s = 1:2
%             % Define Values for Columns
%             model_duct{sim}.([params{p} '_' num2str(s)]) = logninv(P,first_story_columns.([params{p} '_mu_' num2str(s)]),first_story_columns.([params{p} '_beta_' num2str(s)]));
%             duct = model_duct{sim}.([params{p} '_' num2str(s)]);
%             models.([params{p} '_' num2str(s) '_average'])(sim,1) = mean(duct); % standard dev
%             models.([params{p} '_' num2str(s) '_stdev'])(sim,1) = std(duct); % standard dev
%             models.([params{p} '_' num2str(s) '_cov'])(sim,1) = std(duct)/mean(duct); % standard dev
%             models.([params{p} '_' num2str(s) '_beta'])(sim,1) = sqrt(log(std(duct)^2/mean(duct)^2+1)); % lognormal dispersion
%             models.([params{p} '_' num2str(s) '_range'])(sim,1) = (max(duct) - min(duct)); % max difference
%             models.([params{p} '_' num2str(s) '_min'])(sim,1) = min(duct); % max difference
%             models.([params{p} '_' num2str(s) '_max'])(sim,1) = max(duct); % max difference
%             models.([params{p} '_' num2str(s) '_lognorm_range'])(sim,1) = sqrt(log((max(duct)/2 - min(duct)/2)^2/mean(duct)^2+1)); % lognormal normalized max difference
%         end
%     end
% end
% 
% % Save Model Table
% models_table = struct2table(models);
% writetable(models_table, [outputs_dir filesep 'model_table.csv'])
% save([outputs_dir filesep 'model_duct.mat'],'model_duct')

%% Select study models from simulated models
models_table = readtable([outputs_dir filesep 'model_table.csv']);

% Model Ductility Beta Distribution
% histogram(models_table.b_1_cov,'Normalization','pdf')
% m_mean = mean(models_table.b_1_cov);
% m_std = std(models_table.b_1_cov);
% m_mu = log(m_mean^2./sqrt(m_std^2 + m_mean^2));
% m_beta = sqrt(log(m_std^2./m_mean^2+1));
% x_points = 0.01:0.01:1;
% m_pdf = lognpdf(x_points, m_mu, m_beta);
% hold on
% plot(x_points,m_pdf)
% hold off

ranked_data = sortrows(models_table,'b_1_cov');
ranked_data.prob = ((1:1000)/1000)';

% Latin Hypercube Sampling with Targeted Selection
baseline_model_average_duct = mean(first_story_columns.b_hinge_1);
% N = 10;
% for i = 1:N
%     up_prob = i/N;
%     low_prob = up_prob - 1/N;
% %     if i == 10
% %         models_in_range = ranked_data(5*round(ranked_data.b_1_cov/5,2) == 5*round(models2use_betas(i)/5,2),:);
% %     else
% %         models_in_range = ranked_data(round(ranked_data.b_1_cov,2) == round(models2use_betas(i),2),:);
% %     end
%     models_in_range = ranked_data(ranked_data.prob <= up_prob & ranked_data.prob > low_prob,:);
%     [~, idx] = min(abs(models_in_range.b_1_average - baseline_model_average_duct));
%     models2use_id(i) = models_in_range.id(idx);
% end
% models2use = ranked_data(ismember(ranked_data.id,models2use_id),:);
% % Make plots of simulation
% scatter(ranked_data.b_1_cov,ranked_data.prob,'HandleVisibility','off')
% hold on
% scatter(models2use.b_1_cov,models2use.prob,200,'r','filled','p','DisplayName','Selected Models')
% grid on
% xlabel('First-Story Deformation Capacity COV')
% ylabel('P[x>X]')
% legend('Location','northeast')
% fn_format_and_save_plot( outputs_dir, 'Sampled CDF', 4, 1 )
% 
% figure
% scatter(models_table.b_1_cov,models_table.b_1_average,'HandleVisibility','off')
% hold on
% scatter(models2use.b_1_cov,baseline_model_average_duct*ones(N,1),200,'r','filled','p','DisplayName','Selected Models')
% xlabel('First-Story Deformation Capacity COV')
% ylabel('Average First-Story Deformation Capacity')
% grid on
% legend('Location','northeast')
% fn_format_and_save_plot( outputs_dir, 'Sampled Design Space', 4, 1 )
% 
% % Save Selected Model Table
% writetable(models2use, [outputs_dir filesep 'selected_model_table_alt.csv'])

%% Add Uniform and extreme case
models2use = readtable([outputs_dir filesep 'selected_model_table_alt.csv']);
load([outputs_dir filesep 'model_duct.mat'])

% Uniform Case
model_duct{1001}.a_1 = mean(first_story_columns.a_hinge_1);
model_duct{1001}.a_2 = mean(first_story_columns.a_hinge_2);
model_duct{1001}.b_1 = mean(first_story_columns.b_hinge_1);
model_duct{1001}.b_2 = mean(first_story_columns.b_hinge_2);
models2use.id(11) = 1001;
models2use.a_1_beta(11) = 0;
models2use.a_2_beta(11) = 0;
models2use.b_1_beta(11) = 0;
models2use.b_2_beta(11) = 0;

% Extreme Case
highest_cases = ranked_data(end-20:end,:);
[~, idx] = min(abs(highest_cases.b_1_average - baseline_model_average_duct));
model_2_use = highest_cases.id(idx);
models2use(12,:) = ranked_data(ranked_data.id == model_2_use,:);

%% Create plots and tcl script for each model
for m = 11:height(models2use)
    model_dir = [new_models_dir filesep 'Model_' num2str(m)];
    mkdir(model_dir)
    for p = 1:length(params)
        for s = 1:2
            % Write new a and b values to element table
            element.([params{p} '_hinge_' num2str(s)])(ismember(element.id,first_story_columns.id)) = model_duct{models2use.id(m)}.([params{p} '_' num2str(s)]);

            % Make Plots
%             fn_plot_plan_scatter( node_2_use, model_duct{models2use.id(m)}.([params{p} '_' num2str(s)]), model_dir, [params{p} ' value side ' num2str(s)], 0, ['Beta = ' num2str(models2use.([params{p} '_' num2str(s) '_beta'])(m))], [params{p} ' value'], [], 0.1);
        end
    end

    % Copy Model Files
    mkdir([model_dir filesep 'asce_41_data']);
    mkdir([model_dir filesep 'opensees_data']);
    mkdir([model_dir filesep 'model_data']);
    copyfile([os_model_csv_dir filesep 'node.csv'], [model_dir filesep 'model_data'])
    copyfile([os_model_csv_dir filesep 'story.csv'], [model_dir filesep 'model_data'])
    copyfile([os_model_csv_dir filesep 'hinge.csv'], [model_dir filesep 'model_data'])

    % Write Element Analysis Data file with new a and b values
    save([model_dir filesep 'asce_41_data' filesep 'element_analysis.mat'],'element')
    save([model_dir filesep 'asce_41_data' filesep 'joint_analysis.mat'],'joint')

    % Write model.tcl file with new a and b values
    [ ~ ] = fn_define_model( [model_dir filesep 'opensees_data'], node, element, joint, hinge, analysis, '3D', story, [model_dir filesep 'asce_41_data'] );
end

% function [attr_val] = fn_teir_3_model(first_story_columns, attr_name)
%     first_story_columns.idx = ones(height(first_story_columns),1).*(1:24)';
%     teir_1 = [1 6 12 18 19 24];
%     teir_2 = [2 3 4 5 8 14 20 21 22 23];
%     teir_3 = [7 9 10 11 13 15 16 17];
%     attr_val = zeros(height(first_story_columns),1);
%     attr_val(ismember(first_story_columns.idx,teir_1)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_1)));
%     attr_val(ismember(first_story_columns.idx,teir_2)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_2)));
%     attr_val(ismember(first_story_columns.idx,teir_3)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_3)));
% end
