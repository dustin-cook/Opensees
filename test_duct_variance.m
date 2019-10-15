clear all
close all
clc


import asce_41.*
import plotting_tools.*

%% Define Element Properties
ele_prop_table = readtable(['inputs' filesep 'element.csv']);
load(['outputs' filesep 'ICBS_model_3D_fixed' filesep 'NDP_baseline' filesep 'asce_41_data' filesep 'element_analysis.mat'])
load(['outputs' filesep 'ICBS_model_3D_fixed' filesep 'NDP_baseline' filesep 'asce_41_data' filesep 'hinge_analysis.mat'])
load(['outputs' filesep 'ICBS_model_3D_fixed' filesep 'NDP_baseline' filesep 'asce_41_data' filesep 'node_analysis.mat'])

%% Get column hinge properties
first_story_columns = element(element.story == 1 & strcmp(element.type,'column'),:);
for e = 1:height(first_story_columns)
    for s = 1:2
        ele = first_story_columns(e,:);
        ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        ele_side = s;

        % Variance from shear tie spacing
        mean_val = ele_props.(['S_' num2str(ele_side)]);
        standard_deviation = 0.5;
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            ele_props.(['S_' num2str(ele_side)]) = lognrnd(log(mean_val),beta);
            [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, ele_side );
            a.s(i) = hinge.a_hinge;
            b.s(i) = hinge.b_hinge;
        end
        ele_props.(['S_' num2str(ele_side)]) = mean_val;

        % Variance from F'c construction
        mean_val = ele_props.fc_e;
        standard_deviation = 1000;
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            ele_props.fc_e = lognrnd(log(mean_val),beta);
            [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, ele_side );
            a.fc(i) = hinge.a_hinge;
            b.fc(i) = hinge.b_hinge;
        end
        ele_props.fc_e = mean_val; 

        % Variance from F'c Nominal v Expected
        mean_val = mean([ele_props.fc_e,ele_props.fc_n]);
        standard_deviation = ele_props.fc_e - ele_props.fc_n;
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            ele_props.fc_e = lognrnd(log(mean_val),beta);
            [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, ele_side );
            a.fcn(i) = max([hinge.a_hinge,1e-6]);
            b.fcn(i) = max([hinge.b_hinge,1e-6]);
        end
        ele_props.fc_e = mean_val; 

        % Variance from Axial Load
        mean_val = ele.Pmax;
        standard_deviation = abs(ele.Pmax - ele.P_grav);
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            ele.Pmax = lognrnd(log(mean_val),beta);
            [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, ele_side );
            a.axial(i) = hinge.a_hinge;
            b.axial(i) = hinge.b_hinge;
        end
        ele.Pmax = mean_val;

        % Variance from ASCE 41 Testing Data a value
        mean_val = ele.(['a_hinge_' num2str(ele_side)]);
        standard_deviation = abs(mean_val - 0.53*mean_val);
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            a.asce_41(i) = lognrnd(log(mean_val),beta);
        end

        % Variance from ASCE 41 Testing Data b value
        mean_val = ele.(['b_hinge_' num2str(ele_side)]);
        standard_deviation = abs(mean_val - 0.575*mean_val);
        beta = sqrt(log(standard_deviation^2/mean_val^2+1));
        for i = 1:1000
            b.asce_41(i) = lognrnd(log(mean_val),beta);
        end

        % Putting it all together
        mean_s = ele_props.(['S_' num2str(ele_side)]);
        standard_deviation = 0.5;
        beta_s = sqrt(log(standard_deviation^2/mean_s^2+1));

        mean_fc = ele_props.fc_e;
        standard_deviation = 1000;
        beta_fc = sqrt(log(standard_deviation^2/mean_fc^2+1));

        mean_axial = ele.Pmax;
        standard_deviation = abs(ele.Pmax - ele.P_grav);
        beta_axial = sqrt(log(standard_deviation^2/mean_axial^2+1));

        for i = 1:1000
            ele_props.(['S_' num2str(ele_side)]) = lognrnd(log(mean_s),beta_s);
            ele_props.fc_e = lognrnd(log(mean_fc),beta_fc);
            ele.Pmax = lognrnd(log(mean_axial),beta_axial);
            [ hinge, ~, ~ ] = fn_col_hinge( ele, ele_props, ele_side );
            mean_a = hinge.a_hinge;
            mean_b = hinge.b_hinge;

            standard_deviation = abs(mean_a - 0.53*mean_a);
            beta_a = sqrt(log(standard_deviation^2/mean_a^2+1));

            standard_deviation = abs(mean_b - 0.575*mean_b);
            beta_b = sqrt(log(standard_deviation^2/mean_b^2+1));

            a.all(i) = lognrnd(log(mean_a),beta_a);
            b.all(i) = lognrnd(log(mean_b),beta_b);
        end

        %% Fit lognormal curves to each
        clear beta
        terms = {'s','fc','fcn','axial','asce_41','all'};
        for i = 1:length(terms)
            if sum(a.(terms{i})) == 0
                first_story_columns.(['a_mu_' num2str(ele_side) '_' terms{i}])(e,1) = 0;
                first_story_columns.(['a_beta_' num2str(ele_side) '_' terms{i}])(e,1) = 0;
            else
                pHat = lognfit(max(a.(terms{i}),1e-6));
                first_story_columns.(['a_mu_' num2str(ele_side) '_' terms{i}])(e,1) = pHat(1);
                first_story_columns.(['a_beta_' num2str(ele_side) '_' terms{i}])(e,1) = pHat(2);
            end
            if sum(b.(terms{i})) == 0
                first_story_columns.(['b_mu_' num2str(ele_side) '_' terms{i}])(e,1) = 0;
                first_story_columns.(['b_beta_' num2str(ele_side) '_' terms{i}])(e,1) = 0;
            else
                pHat = lognfit(max(b.(terms{i}),1e-6));
                first_story_columns.(['b_mu_' num2str(ele_side) '_' terms{i}])(e,1) = pHat(1);
                first_story_columns.(['b_beta_' num2str(ele_side) '_' terms{i}])(e,1) = pHat(2);
            end
        end
    end
end

% Create plots to illustrate to abbie
% Plan view plots
node_2_use = node(ismember(node.id,first_story_columns.node_1),:);
fn_plot_plan_scatter( node_2_use, first_story_columns.b_beta_1_all, pwd, 'temp_plot', 0, 'First Story Column Base', 'Total Lognormal Dispersion', [], max(first_story_columns.b_beta_1_all));

% Single element PDF plots
