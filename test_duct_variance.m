clear all
close all
clc

import asce_41.*
import plotting_tools.*
import opensees.write_tcl.*

% Define Model
analysis.model_id = 12;
analysis.proceedure = 'NDP';
analysis.id = 'ida_new';

% Secondary options
analysis.type = 1;
analysis.nonlinear = 1;
analysis.stories_nonlinear = 3; % Default to all modeling all stories as nonlinear when doing NDP
analysis.model_type = 2; % 1 = SDOF, 2 = MDOF (default)
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;
analysis.suppress_outputs = 1;
analysis.joint_explicit = 1;
analysis.joint_model = 1; % 1 = beam/column elements, 2 = joint 3D

%% Define inputs and outputs directories
analysis_name = [analysis.proceedure '_' analysis.id];
inputs_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep analysis_name filesep 'asce_41_data'];
os_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep analysis_name filesep 'opensees_data'];
os_model_csv_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep analysis_name filesep 'model_data'];
pushover_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep analysis_name filesep 'pushover'];
outputs_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep analysis_name filesep 'sensitivity_study'];
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

% Overwrite Pmax from Dynamic with Pmax from Pushover
pushover = load([pushover_dir filesep 'element_analysis.mat']);
element.Pmax = pushover.element.Pmax;

% Filter to just first story columns
first_story_columns = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,first_story_columns.node_1),:);

% Set parameters to investigate
params = {'a', 'b'};
model_name = {'fc', 's', 'axial', 'asce_41','all'};
model_descr = {'Concrete Strength','Shear Tie Spacing','Axial Load','ASCE 41 Uncertainty','Cumulative Uncertianty'};
model_type = [1 1 1 2 3];

%% Get column hinge properties
% for e = 1:height(first_story_columns)
%     e
%     for s = 1:2
%         % Define Model Parameters
%         ele = first_story_columns(e,:);
%         ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
%         model_variable_1 = {'ele_props.fc_e', ['ele_props.S_' num2str(s)], 'ele.Pmax', ['ele.a_hinge_' num2str(s)]};
%         model_variable_2 = {'ele_props.fc_e', ['ele_props.S_' num2str(s)], 'ele.Pmax', ['ele.b_hinge_' num2str(s)]};
%         model_mean_1 = [ele_props.fc_e, ele_props.(['S_' num2str(s)]), ele.Pmax, ele.(['a_hinge_' num2str(s)])];
%         model_mean_2 = [ele_props.fc_e, ele_props.(['S_' num2str(s)]), ele.Pmax, ele.(['b_hinge_' num2str(s)])];
%         model_std_1 = [1000, 0.5, abs(ele.Pmax - ele.P_grav), ele.(['a_hinge_' num2str(s)])*(1 - 0.53)];
%         model_std_2 = [1000, 0.5, abs(ele.Pmax - ele.P_grav), ele.(['b_hinge_' num2str(s)])*(1 - 0.575)];
%         model_mu_1 = log(model_mean_1.^2./sqrt(model_std_1.^2 + model_mean_1.^2));
%         model_mu_2 = log(model_mean_2.^2./sqrt(model_std_2.^2 + model_mean_2.^2));
%         model_beta_1 = sqrt(log(model_std_1.^2./model_mean_1.^2+1));
%         model_beta_2 = sqrt(log(model_std_2.^2./model_mean_2.^2+1));
%         
%         % Sample Hinge Calculation with model uncertainty
%         for m = 1:length(model_name)
%             % Element Properties
%             ele = first_story_columns(e,:);
%             ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
%             
%             % Sample Values
%             if model_type(m) == 1
%                 for i = 1:1000
%                     eval(sprintf('%s = lognrnd(model_mu_1(m),model_beta_1(m));',model_variable_1{m}))
%                     [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, s );
%                     samples.a.([model_name{m} '_' num2str(s)])(e,i) = hinge.a_hinge;
%                     samples.b.([model_name{m} '_' num2str(s)])(e,i) = hinge.b_hinge;
%                 end
%             elseif model_type(m) == 2
%                 for i = 1:1000
%                     if model_mean_1(m) == 0
%                         samples.a.([model_name{m} '_' num2str(s)])(e,i) = 0;
%                     else
%                         samples.a.([model_name{m} '_' num2str(s)])(e,i) = lognrnd(model_mu_1(m),model_beta_1(m));
%                     end
%                     if model_mean_2(m) == 0
%                         samples.b.([model_name{m} '_' num2str(s)])(e,i) = 0;
%                     else
%                         samples.b.([model_name{m} '_' num2str(s)])(e,i) = lognrnd(model_mu_2(m),model_beta_2(m));
%                     end
%                 end
%             elseif model_type(m) == 3
%                 for i = 1:1000
%                     for j = 1:length(model_variable_1)-1 % take all but last entry which is the asce 41 entry
%                         eval(sprintf('%s = lognrnd(model_mu_1(j),model_beta_1(j));',model_variable_1{j}))
%                     end
%                     [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, s );
%                     a_std = hinge.a_hinge*(1 - 0.53);
%                     b_std = hinge.b_hinge*(1 - 0.575);
%                     a_mu = log(hinge.a_hinge^2./sqrt(a_std^2 + hinge.a_hinge^2));
%                     b_mu = log(hinge.b_hinge^2./sqrt(b_std^2 + hinge.b_hinge^2));
%                     a_beta = max([sqrt(log(a_std^2/hinge.a_hinge^2+1)),0]);
%                     b_beta = max([sqrt(log(b_std^2/hinge.b_hinge^2+1)),0]);
%                     if hinge.a_hinge == 0
%                         samples.a.([model_name{m} '_' num2str(s)])(e,i) = 0;
%                     else
%                         samples.a.([model_name{m} '_' num2str(s)])(e,i) = lognrnd(a_mu,a_beta);
%                     end
%                     if hinge.b_hinge == 0
%                         samples.b.([model_name{m} '_' num2str(s)])(e,i) = 0;
%                     else
%                         samples.b.([model_name{m} '_' num2str(s)])(e,i) = lognrnd(b_mu,b_beta);
%                     end
%                 end
%             end
%             
%             for p = 1:length(params)
%                 % Sort Samples
%                 samples.(params{p}).([model_name{m} '_' num2str(s)])(e,:) = sort(samples.(params{p}).([model_name{m} '_' num2str(s)])(e,:),2);
% 
%                 % Fit lognormal curves to each
%                 mean_val = mean(samples.(params{p}).([model_name{m} '_' num2str(s)])(e,:));
%                 std_val = std(samples.(params{p}).([model_name{m} '_' num2str(s)])(e,:));
%                 first_story_columns.([params{p} '_mean_' num2str(s) '_' model_name{m}])(e,1) = mean_val;
%                 first_story_columns.([params{p} '_std_' num2str(s) '_' model_name{m}])(e,1) = std_val;
%                 if mean_val == 0
%                     first_story_columns.([params{p} '_mu_' num2str(s) '_' model_name{m}])(e,1) = 0;
%                     first_story_columns.([params{p} '_beta_' num2str(s) '_' model_name{m}])(e,1) = 0;
%                 else
%                     first_story_columns.([params{p} '_mu_' num2str(s) '_' model_name{m}])(e,1) = log(mean_val^2/sqrt(std_val^2 + mean_val^2));
%                     first_story_columns.([params{p} '_beta_' num2str(s) '_' model_name{m}])(e,1) = sqrt(log(std_val^2/mean_val^2+1));
%                 end
%             end
%         end
%     end
% end
% 
% % Save Simulated Data
% save([outputs_dir filesep 'sim_data.mat'],'samples','first_story_columns')

%% Post Process Samplings
load([outputs_dir filesep 'sim_data.mat'])

% % Create plots of pds for each element
% x_points = 0.0001:0.0001:0.05;
% pdf_plot_dir = [outputs_dir filesep 'element_pdf_plots'];
% mkdir(pdf_plot_dir)
% for e = 1:height(first_story_columns)
%     for s = 1:2
%         for p = 1:length(params)
%             hold on
%             for m = 1:length(model_name)
%                 pdf = lognpdf(x_points,first_story_columns.([params{p} '_mu_' num2str(s) '_' model_name{m}])(e),first_story_columns.([params{p} '_beta_' num2str(s) '_' model_name{m}])(e));
%                 plot(x_points,pdf./max(pdf),'DisplayName',model_descr{m})
%             end
%             xlabel([params{p} ' hinge value (rot)'])
%             ylabel('Normalized Probaility Density')
%             grid on
%             fn_format_and_save_plot( pdf_plot_dir, [params{p} '_' num2str(s) '_ele_' num2str(e)], 3, 1 )
%         end
%     end
% end

% Define Models Run
model_name = {'Baseline', 'Uniform', '3 Tier', 'fc', 's', 'axial', 'asce_41','all'};
model_sims = [1 1 1 3 3 3 3 3];
model_var = {'first_story_columns.([params{p} ''_hinge_'' num2str(s)])', ...
             'mean(first_story_columns.([params{p} ''_hinge_'' num2str(s)]))*ones(height(first_story_columns),1)',...
             'fn_teir_3_model(first_story_columns, [params{p} ''_hinge_'' num2str(s)])',...
             'diag(samples.(params{p}).([model_name{m} ''_'' num2str(s)])(1:height(first_story_columns),idx))',...
             'diag(samples.(params{p}).([model_name{m} ''_'' num2str(s)])(1:height(first_story_columns),idx))',...
             'diag(samples.(params{p}).([model_name{m} ''_'' num2str(s)])(1:height(first_story_columns),idx))',...
             'diag(samples.(params{p}).([model_name{m} ''_'' num2str(s)])(1:height(first_story_columns),idx))',...
             'diag(samples.(params{p}).([model_name{m} ''_'' num2str(s)])(1:height(first_story_columns),idx))'};

% Go through each model, calc diversity, write tcl files, and plot results
m_id = 0;
for m = 1:length(model_name)
    model_dir = [outputs_dir filesep model_name{m}];
    mkdir(model_dir)
    for sim = 1:model_sims(m)
        m_id = m_id + 1;
        models.id(m_id,1) = m_id;
        models.type{m_id,1} = model_name{m};
        idx = randi(1000,height(first_story_columns));
        for p = 1:length(params)
            for s = 1:2
                % Define Values for Columns
                duct_val = eval(model_var{m});

                % Measure Ductility Parameters
                models.([params{p} '_' num2str(s) '_average'])(m_id,1) = mean(duct_val); % standard dev
                models.([params{p} '_' num2str(s) '_diversity_1'])(m_id,1) = std(duct_val); % standard dev
                models.([params{p} '_' num2str(s) '_diversity_2'])(m_id,1) = sqrt(log(std(duct_val)^2/mean(duct_val)^2+1)); % lognormal dispersion
                models.([params{p} '_' num2str(s) '_diversity_3'])(m_id,1) = (max(duct_val) - min(duct_val))/2; % max difference
                models.([params{p} '_' num2str(s) '_diversity_4'])(m_id,1) = sqrt(log((max(duct_val)/2 - min(duct_val)/2)^2/mean(duct_val)^2+1)); % lognormal normalized max difference
                
                % Write new a and b values to element table
                element.([params{p} '_hinge_' num2str(s)])(ismember(element.id,first_story_columns.id)) = duct_val;
                
                % Make Plots
%                 fn_plot_plan_scatter( node_2_use, duct_val, model_dir, [params{p} ' value side ' num2str(s) ' ' model_name{m} ' sim ' num2str(sim)], 0, ['Element Lognormal Dispersion = ' num2str(models.([params{p} '_' num2str(s) '_diversity_2'])(m_id,1))], [params{p} ' value'], [], max(duct_val));
            end
        end
        % Copy Model Files
        model_dir = [new_models_dir filesep model_name{m} '_' num2str(sim)];
        mkdir(model_dir);
        mkdir([model_dir filesep 'asce_41_data']);
        mkdir([model_dir filesep 'opensees_data']);
        mkdir([model_dir filesep 'model_data']);
        copyfile([os_model_csv_dir filesep 'node.csv'], [model_dir filesep 'model_data'])
        copyfile([os_model_csv_dir filesep 'story.csv'], [model_dir filesep 'model_data'])
        copyfile([os_model_csv_dir filesep 'hinge.csv'], [model_dir filesep 'model_data'])
                
        % Write Element Analysis Data file with new a and b values
        save([model_dir filesep 'asce_41_data' filesep 'element_analysis.mat'],'samples','element')
        save([model_dir filesep 'asce_41_data' filesep 'joint_analysis.mat'],'samples','joint')

        % Write model.tcl file with new a and b values
        [ ~ ] = fn_define_model( [model_dir filesep 'opensees_data'], node, element, joint, hinge, analysis, '3D', story, [model_dir filesep 'asce_41_data'] );
    end
end

% Save Model Table
writetable(struct2table(models), [outputs_dir filesep 'model_table.csv'])


function [attr_val] = fn_teir_3_model(first_story_columns, attr_name)
    first_story_columns.idx = ones(height(first_story_columns),1).*(1:24)';
    teir_1 = [1 6 12 18 19 24];
    teir_2 = [2 3 4 5 8 14 20 21 22 23];
    teir_3 = [7 9 10 11 13 15 16 17];
    attr_val = zeros(height(first_story_columns),1);
    attr_val(ismember(first_story_columns.idx,teir_1)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_1)));
    attr_val(ismember(first_story_columns.idx,teir_2)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_2)));
    attr_val(ismember(first_story_columns.idx,teir_3)) = mean(first_story_columns.(attr_name)(ismember(first_story_columns.idx,teir_3)));
end
