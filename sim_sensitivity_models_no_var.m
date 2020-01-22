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
model_name = {'ductility_no_var'};
params{1} = {'a_hinge', 'b_hinge'};
dev_params{1} = [0.5264, 0.5920];
prime_param = {'b_hinge_1'};
max_val = [0.1];

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
    
    %% Create array of models from low to high parameter range
    for i = 1:6 % six new models
        for e = 1:height(first_story_columns)
            for p = 1:length(params{m})
                for s = 1:2
                    i_percentile = normcdf(i-3.5);
                    ele_mu = first_story_columns.([params{m}{p} '_mu_' num2str(s)])(e,:);
                    ele_beta = first_story_columns.([params{m}{p} '_beta_' num2str(s)])(e,:);
                    model_param_vals{i}.([params{m}{p} '_' num2str(s)])(e,:) = logninv(i_percentile,ele_mu,ele_beta);
                end
            end
        end
    end
    
    %% Create plots and tcl script for each model
    for mod = 1:length(model_param_vals)
        model_dir = [new_models_dir filesep 'model_' num2str(mod)];
        mkdir(model_dir)
        for p = 1:length(params{m})
            for s = 1:2
                % Write new parameter values to element table
                element.([params{m}{p} '_' num2str(s)])(ismember(element.id,first_story_columns.id)) = model_param_vals{mod}.([params{m}{p} '_' num2str(s)]);

                % Make Plots
                fn_plot_plan_scatter( node_2_use, model_param_vals{mod}.([params{m}{p} '_' num2str(s)]), model_dir, [params{m}{p} ' value side ' num2str(s)], 0, ['Mean = ' num2str(mean(model_param_vals{mod}.([params{m}{p} '_' num2str(s)])))], [params{m}{p}], [], max_val(m));
            end
        end
        
        % Write CP values based on 0.7 b from table 10-9
        element.cp_1(ismember(element.id,first_story_columns.id)) = 0.7*element.b_hinge_1(ismember(element.id,first_story_columns.id));
        element.cp_2(ismember(element.id,first_story_columns.id)) = 0.7*element.b_hinge_2(ismember(element.id,first_story_columns.id));

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



