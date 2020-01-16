clear all
close all
clc

% Creates variations on pre-defined models 
% Only for 2D models
% Only modifies the first story columns

import asce_41.*
import plotting_tools.*
import opensees.write_tcl.*

rng(111389)

% Define Model
analysis.model_id = 18;
analysis.proceedure = 'NDP';
analysis.id = 'baseline';
analysis.model = 'ICBS_model_5ew_col_base';

% Secondary options
analysis.type = 1;
analysis.nonlinear = 1;
analysis.stories_nonlinear = 1; % Default to all modeling all stories as nonlinear when doing NDP
analysis.elastic_beams = 1;
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
    
%% Define sensitivity study parameters
model_name = {'ductility', 'strength'};
params{1} = {'a_hinge', 'b_hinge'};
params{2} = {'Mn_pos', 'Mn_neg', 'Mp_pos', 'Mp_neg'};
dev_params{1} = [0.5264, 0.5920];
dev_params{2} = [0.75, 0.75, 0.75, 0.75];
prime_param = {'b_hinge_1', 'Mn_pos_1'};
max_val = [0.1, 2e7];

%% Create Models for each sensitivity study
for m = 1:length(model_name)
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

    % Create outputs dir
    outputs_dir = ['outputs' filesep analysis.model filesep analysis_name filesep 'sensitivity_study' filesep model_name{m}];
    new_models_dir = [outputs_dir filesep 'model_files'];
    mkdir(outputs_dir)
    mkdir(new_models_dir)

    %% Get column hinge properties
    for e = 1:height(first_story_columns)
        for s = 1:2
            for p = 1:length(params{m})
                % Define Model Parameters
                ele = first_story_columns(e,:);
                ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
                first_story_columns.([params{m}{p} '_mean_' num2str(s)])(e,:) = ele.([params{m}{p} '_' num2str(s)]);
                first_story_columns.([params{m}{p} '_std_' num2str(s)])(e,:) = ele.([params{m}{p} '_' num2str(s)])*(1 - dev_params{m}(p));
                first_story_columns.([params{m}{p} '_mu_' num2str(s)])(e,:)= log(ele.([params{m}{p} '_' num2str(s)]).^2./sqrt(first_story_columns.([params{m}{p} '_std_' num2str(s)])(e,:).^2 + ele.([params{m}{p} '_' num2str(s)]).^2));
                first_story_columns.([params{m}{p} '_beta_' num2str(s)])(e,:) = sqrt(log(first_story_columns.([params{m}{p} '_std_' num2str(s)])(e,:).^2./ele.([params{m}{p} '_' num2str(s)]).^2+1));
            end
        end
    end

    % Save Simulated Data
    save([outputs_dir filesep 'sim_data.mat'],'first_story_columns')

    %% Simulate 1000 models
    for sim = 1:1000
        models{m}.id(sim,1) = sim;
        P = rand(height(first_story_columns),1);
        for p = 1:length(params{m})
            for s = 1:2
                % Define Values for Columns
                model_param_vals{m}{sim}.([params{m}{p} '_' num2str(s)]) = logninv(P,first_story_columns.([params{m}{p} '_mu_' num2str(s)]),first_story_columns.([params{m}{p} '_beta_' num2str(s)]));
                param_val = model_param_vals{m}{sim}.([params{m}{p} '_' num2str(s)]);
                models{m}.([params{m}{p} '_' num2str(s) '_average'])(sim,1) = mean(param_val); % mean
                models{m}.([params{m}{p} '_' num2str(s) '_stdev'])(sim,1) = std(param_val); % standard dev
                models{m}.([params{m}{p} '_' num2str(s) '_cov'])(sim,1) = std(param_val)/mean(param_val); % coefficient of variatio
                models{m}.([params{m}{p} '_' num2str(s) '_beta'])(sim,1) = sqrt(log(std(param_val)^2/mean(param_val)^2+1)); % lognormal dispersion
                models{m}.([params{m}{p} '_' num2str(s) '_range'])(sim,1) = (max(param_val) - min(param_val)); % max difference
                models{m}.([params{m}{p} '_' num2str(s) '_min'])(sim,1) = min(param_val); % max difference
                models{m}.([params{m}{p} '_' num2str(s) '_max'])(sim,1) = max(param_val); % max difference
                models{m}.([params{m}{p} '_' num2str(s) '_lognorm_range'])(sim,1) = sqrt(log((max(param_val)/2 - min(param_val)/2)^2/mean(param_val)^2+1)); % lognormal normalized max difference
            end
        end
    end
    
    % Save Model Table
    models_table = struct2table(models{m});
    writetable(models_table, [outputs_dir filesep 'model_table.csv'])
    save([outputs_dir filesep 'param_vals.mat'],'param_val')
    
    %% Select study models from simulated models
    models_table = readtable([outputs_dir filesep 'model_table.csv']);
    
    % Model Ductility Beta Distribution
    histogram(models_table.([prime_param{m} '_cov']),'Normalization','pdf')
    m_mean = mean(models_table.([prime_param{m} '_cov']));
    m_std = std(models_table.([prime_param{m} '_cov']));
    m_mu = log(m_mean^2./sqrt(m_std^2 + m_mean^2));
    m_beta = sqrt(log(m_std^2./m_mean^2+1));
    x_points = 0.01:0.01:1;
    m_pdf = lognpdf(x_points, m_mu, m_beta);
    hold on
    plot(x_points,m_pdf)
    hold off
    fn_format_and_save_plot( outputs_dir, 'Simulated PDF', 4, 1 )
    
    % Latin Hypercube Sampling with Targeted Selection
    ranked_data = sortrows(models_table,[prime_param{m} '_cov']);
    ranked_data.prob = ((1:1000)/1000)';
    baseline_model_average_param = mean(first_story_columns.([prime_param{m}]));
    baseline_model_cov_param = std(first_story_columns.([prime_param{m}]))/mean(first_story_columns.([prime_param{m}]));
    N = 5;
    for i = 1:N
        up_prob = i/N;
        low_prob = up_prob - 1/N;
        models_in_range = ranked_data(ranked_data.prob <= up_prob & ranked_data.prob > low_prob,:);
        [~, idx] = min(abs(models_in_range.([prime_param{m} '_average']) - baseline_model_average_param));
        models2use_id(i) = models_in_range.id(idx);
    end
    models2use = ranked_data(ismember(ranked_data.id,models2use_id),:);
    
    % Make plots of simulation
    scatter(ranked_data.([prime_param{m} '_cov']),ranked_data.prob,'HandleVisibility','off')
    hold on
    scatter(models2use.([prime_param{m} '_cov']),models2use.prob,200,'r','filled','p','DisplayName','Selected Models')
    grid on
    xlabel(['First-Story ' prime_param{m} ' COV'])
    ylabel('P[x>X]')
    legend('Location','northeast')
    fn_format_and_save_plot( outputs_dir, 'Sampled CDF', 4, 1 )

    figure
    scatter(models_table.([prime_param{m} '_cov']),models_table.([prime_param{m} '_average']),'HandleVisibility','off')
    hold on
    scatter(models2use.([prime_param{m} '_cov']),baseline_model_average_param*ones(N,1),200,'r','filled','p','DisplayName','Selected Models')
    xlabel(['First-Story ' prime_param{m} ' COV'])
    ylabel(['Average First-Story' prime_param{m}])
    grid on
    legend('Location','northeast')
    fn_format_and_save_plot( outputs_dir, 'Sampled Design Space', 4, 1 )

    % Save Selected Model Table
    writetable(models2use, [outputs_dir filesep 'selected_model_table_alt.csv'])
    
    %% Add Uniform case
    models2use = readtable([outputs_dir filesep 'selected_model_table_alt.csv']);
    load([outputs_dir filesep 'param_vals.mat'])
    models2use.id(6) = 1001;
    for p = 1:length(params{m})
        for s = 1:2
            model_param_vals{m}{1001}.([params{m}{p} '_' num2str(s)]) = mean(first_story_columns.([params{m}{p} '_' num2str(s)]))*ones(length(first_story_columns.([params{m}{p} '_' num2str(s)])),1);
            models2use.([params{m}{p} '_' num2str(s) '_beta'])(6) = 0;
        end
    end

    % Extreme Case
    highest_cases = ranked_data(end-20:end,:);
    [~, idx] = min(abs(highest_cases.([prime_param{m} '_average']) - baseline_model_average_param));
    model_2_use = highest_cases.id(idx);
    models2use(7,:) = ranked_data(ranked_data.id == model_2_use,:);
    
    %% Create plots and tcl script for each model
    for mod = 1:height(models2use)
        model_dir = [new_models_dir filesep 'model_' num2str(mod)];
        mkdir(model_dir)
        for p = 1:length(params{m})
            for s = 1:2
                % Write new a and b values to element table
                element.([params{m}{p} '_' num2str(s)])(ismember(element.id,first_story_columns.id)) = model_param_vals{m}{models2use.id(mod)}.([params{m}{p} '_' num2str(s)]);

                % Make Plots
                fn_plot_plan_scatter( node_2_use, model_param_vals{m}{models2use.id(mod)}.([params{m}{p} '_' num2str(s)]), model_dir, [params{m}{p} ' value side ' num2str(s)], 0, ['Beta = ' num2str(models2use.([params{m}{p} '_' num2str(s) '_beta'])(mod))], [params{m}{p} ' value'], [], max_val(m));
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

        % Write model.tcl file with new values
        [ ~ ] = fn_define_model( [model_dir filesep 'opensees_data'], node, element, joint, hinge, analysis, '2D', story, [model_dir filesep 'asce_41_data'] );
    end
end



