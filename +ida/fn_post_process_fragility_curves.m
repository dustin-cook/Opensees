function [ ] = fn_post_process_fragility_curves(analysis,model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Define Write Directory
write_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Compare'];
mkdir(write_dir)

% Load Fragility Data
read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
load([read_dir filesep 'frag_curves.mat'])
load([read_dir filesep 'stripe_table.mat'])
ida_table = readtable([read_dir filesep 'ida_table.csv']);
ida_table_full = readtable([read_dir filesep 'ida_table_full.csv']);

% Define mechanisms and metrics to compare
mechanism = {'cols_walls_1','cols_walls_1','drift','drift','drift','drift','drift','gravity','gravity','gravity','gravity','gravity','adjacent_comp','adjacent_comp','adjacent_comp'};
params = {'cp','b_e','idr_1','idr_2','idr_3','idr_4','idr_5','percent_lost_10','percent_lost_25','percent_lost_50','percent_lost_75','percent_lost_100','any','any_frame','all'};

% Use Collapse as consequence metrics
consequence_theta = frag_curves.collapse.theta;
consequence_beta = frag_curves.collapse.beta;
collapse_target = ida_table_full.collapse;

%% Do some tests
% x_points = 0:0.01:2;
% hold on
% plot(x_points,logncdf(x_points,log(consequence_theta),consequence_beta))
% scatter(stripe.asce41_sa_ew,stripe.num_collapse_tot ./ stripe.num_gms)
% gms = unique(ida_table.gm_name);
% idx = 1;
% for gm = 1:length(gms)
%     gm_collapse_sas = ida_table.sa_x(ida_table.collapse > 0 & strcmp(ida_table.gm_name,gms{gm}));
%     if ~isempty(gm_collapse_sas)
%         collapse_sas(idx) = min(gm_collapse_sas);
%         idx = idx + 1;
%     end
% end
% collapse_sa_dist = sort(collapse_sas);
% collapse_sa_rank = (1:length(collapse_sa_dist)) / length(collapse_sa_dist);
% plot(collapse_sa_dist,collapse_sa_rank)
% scatter(ida_table_full.sa_x,collapse_target)
% B = mnrfit(log(ida_table_full.sa_x),categorical(~collapse_target));
% % B = glmfit(log(ida_table_full.sa_x),collapse_target,'binomial','link','probit');
% z_val = B(1) + log(x_points)*B(2);
% prob_val = 1 ./ (1+exp(-z_val));
% plot(x_points,prob_val)
% 
% prediction = zeros(length(collapse_target),1);
% prediction(prob_val > 0.5) = 1;
% LR_accuracy = sum(collapse_target == prediction)/length(collapse_target);
% TP = sum(and(prediction,collapse_target));
% TN = sum(and(~prediction,~collapse_target));
% FP = sum(and(prediction,~collapse_target));
% FN = sum(and(~prediction,collapse_target));
% LR_precision = TP / (TP+FP);
% LR_recall = TP / (FN+TP);
% LR_F1 = 2*LR_precision*LR_recall / (LR_precision + LR_recall);

%% For each component metrics, quantify comparison
% metrics = table;
% idx = 0;
% for p = 1:length(params)
%     for n = 1:length(frag_curves.(mechanism{p}).(params{p}).theta)
%         idx = idx + 1;
%         metrics.id(idx,1) = idx;
%         metrics.mechanism{idx,1} = mechanism{p};
%         metrics.parameter{idx,1} = params{p};
%         if ~strcmp(mechanism{p},'cols_walls_1')
%             metrics.num_components(idx,1) = NaN;
%             metrics.percent_mechanism(idx,1) = NaN;
%             plot_dir = [write_dir filesep mechanism{p} '_' params{p}];
%         else
%             metrics.num_components(idx,1) = frag_curves.(mechanism{p}).(params{p}).num_comp(n);
%             metrics.percent_mechanism(idx,1) = frag_curves.(mechanism{p}).(params{p}).prct_mech(n);
%             plot_dir = [write_dir filesep mechanism{p} '_' params{p} '_' num2str(metrics.num_components(idx,1)) '_comps'];
%         end
%         component_theta = frag_curves.(mechanism{p}).(params{p}).theta(n);
%         component_beta = frag_curves.(mechanism{p}).(params{p}).beta(n);
%         [metrics_data(idx,:)] = fn_comparison_metrics(consequence_theta, consequence_beta, component_theta, component_beta, plot_dir);
%     end
% end
% 
% % Write Table of Metrics
% metrics_table = [metrics, metrics_data];
% writetable(metrics_table,[write_dir filesep 'metrics_table.csv'])

%% For each component metric, quantify comparison logistic regression
mechanism = {'cols_walls_1','cols_walls_1','drift', 'gravity'};
params = {'cp','b_e','idr','percent_lost'};
LR_metrics = table;
idx = 0;
for p = 1:length(params)
    if strcmp(mechanism{p},'drift')
        LR_parameters = {'drift'};
        LR_x_names = {'Max Story Drift'};
    elseif strcmp(mechanism{p},'gravity')
        LR_parameters = {'gravity'};
        LR_x_names = {'Ratio of Gravity Capacity Lost'};
    else
        LR_parameters = {'max', 'mean', 'num', 'percent'};
        LR_x_names = {'Maximum Component Ratio of Demand to CP', 'Average Component Ratio of Demand to CP', 'Number of Components Exceeding CP', 'Percent of Mechanism Exceeding CP'};
    end
    for lr = 1:length(LR_parameters)
        idx = idx + 1;
        LR_metrics.id(idx,1) = idx;
        LR_metrics.mechanism{idx,1} = mechanism{p};
        LR_metrics.parameter{idx,1} = params{p};
        LR_metrics.metric{idx,1} = LR_parameters{lr};
        if strcmp(mechanism{p},'drift')
            component_response_param = max(ida_table_full.drift_x,ida_table_full.drift_z);
            plot_dir = [write_dir filesep 'LR_' mechanism{p}];
        elseif strcmp(mechanism{p},'gravity')
            component_response_param = ida_table_full.gravity_load_lost_ratio;
            plot_dir = [write_dir filesep 'LR_' mechanism{p}];
        else
            component_response_param = ida_table_full.([mechanism{p} '_' LR_parameters{lr} '_' params{p}]);
            plot_dir = [write_dir filesep 'LR_' mechanism{p} '_' params{p} '_' LR_parameters{lr}];
        end
        [LR_metrics_data(idx,:)] = fn_logistic_regression(component_response_param, collapse_target, plot_dir, LR_x_names{lr});
    end
end

% Write Table of Metrics
LR_metrics_table = [LR_metrics, LR_metrics_data];
writetable(LR_metrics_table,[write_dir filesep 'LR_metrics_table.csv'])

end

function [metrics] = fn_comparison_metrics(consequence_theta, consequence_beta, component_theta, component_beta, plot_dir)
import plotting_tools.fn_format_and_save_plot

%% Define Consequence Curve
x_points = 0:0.01:4;
consequence_pdf = lognpdf(x_points,log(consequence_theta),consequence_beta);
component_pdf = lognpdf(x_points,log(component_theta),component_beta);
consequence_cdf = logncdf(x_points,log(consequence_theta),consequence_beta);
component_cdf = logncdf(x_points,log(component_theta),component_beta);

%% Reliability Index
[Mu_C,V_C] = lognstat(log(consequence_theta),consequence_beta);
[Mu_D,V_D] = lognstat(log(component_theta),component_beta);
metrics.RI = (Mu_C - Mu_D)/sqrt(V_C + V_D);

%% Probability of Failure
metrics.POF = normcdf(-metrics.RI);
hold on
plot(x_points,component_pdf,'DisplayName','ASCE 41 CP')
plot(x_points,consequence_pdf,'DisplayName','Collapse')
xlabel('Sa (T_{EW}) (g)')
ylabel('Probability Density')
fn_format_and_save_plot( plot_dir, 'Probabilty of Failure', 3, 1 )

%% Ratio of Median Values
metrics.CMR = consequence_theta ./ component_theta;

%% Difference of Median Values
metrics.abs_diff_sa = consequence_theta - component_theta;
metrics.percent_diff = (consequence_theta - component_theta) / component_theta;

%% Difference of Percentiles (90th to 10th)
[~, idx_10th] = min(abs(consequence_cdf - 0.1));
[~, idx_90th] = min(abs(component_cdf - 0.9));
metrics.abs_diff_90_10 = x_points(idx_10th) - x_points(idx_90th);

%% P[Consequence|Median Component Response] (maybe plot for a whole range of percentiles)
probs = 0.01:0.01:0.99;
for p = 1:length(probs)
    [~, idx_comp_prob] = min(abs(component_cdf - probs(p)));
    prob_consequece_given_component_percentile(p) = consequence_cdf(idx_comp_prob);
end
metrics.prob_consequece_given_component_med = prob_consequece_given_component_percentile(probs == 0.5);
plot(probs,prob_consequece_given_component_percentile)
xlabel('Probability of Component Exceeding CP')
ylabel('Probability of Collapse')
fn_format_and_save_plot( plot_dir, 'Probability QQ', 3, 1 )

% Change Metrics to table for outputs
metrics = struct2table(metrics);

end

function [metrics] = fn_logistic_regression(max_component_response, collapse_target, plot_dir, x_label)
import plotting_tools.fn_format_and_save_plot

%% Logistic Regression
% ida_table = readtable([read_dir filesep 'ida_table.csv']);
% collapse_target = ida_table.collapse;
% max_component_response = ida_table.max_cp_cols_walls_1;
B = mnrfit(log(max_component_response),categorical(~collapse_target));
z_val = B(1) + log(max_component_response)*B(2);
prob_val = 1 ./ (1+exp(-z_val));
metrics.LR_value_coef = B(1);
metrics.LR_correlation_coef = B(2);
metrics.LR_threshold_response_ratio = -B(1)/B(2);
metrics.LR_prob_consequece_given_component = 1 / (1+exp(-(B(1) + 1*B(2))));
prediction = zeros(length(collapse_target),1);
prediction(prob_val > 0.5) = 1;
metrics.LR_accuracy = sum(collapse_target == prediction)/length(collapse_target);
TP = sum(and(prediction,collapse_target));
TN = sum(and(~prediction,~collapse_target));
FP = sum(and(prediction,~collapse_target));
FN = sum(and(~prediction,collapse_target));
metrics.LR_precision = TP / (TP+FP);
metrics.LR_recall = TP / (FN+TP);
metrics.LR_F1 = 2*metrics.LR_precision*metrics.LR_recall / (metrics.LR_precision + metrics.LR_recall);

hold on
scatter(max_component_response,collapse_target,'DisplayName','Recorded Collapse')
scatter(max_component_response,prob_val,'DisplayName','Predicted Collapse')
xlabel(x_label)
ylabel('Probability of Collapse')
fn_format_and_save_plot( plot_dir, 'Logistic Regression', 3, 1 )

% Change Metrics to table for outputs
metrics = struct2table(metrics);

end
