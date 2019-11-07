function [ ] = fn_LR_classification(analysis,model,gm_set_table)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Load Fragility Data
read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
ida_table = readtable([read_dir filesep 'ida_table.csv']);

%% Fake new collapse cases above for each ground motion such that LR is balanced (or just take the just before collapse and just after collapse for each)
% for gm = 1:height(gm_set_table)
%     gm_response_no_col = ida_table(strcmp(ida_table.eq_name,gm_set_table.eq_name(gm)) & ida_table.collapse == 0,:);
%     Sa_jbc = max(gm_response_no_col.sa_x);
%     ida_table.jbc(strcmp(ida_table.eq_name,gm_set_table.eq_name(gm)) & ida_table.sa_x == Sa_jbc) = 1;
% end
% filt_ida_table = ida_table(ida_table.jbc == 1 | ida_table.collapse > 0,:);

new_ida_table = table;
for gm = 1:height(gm_set_table)
    gm_response = ida_table(strcmp(ida_table.eq_name,gm_set_table.eq_name(gm)),:);
    gm_response.collapse(end) = 1;
    no_col_cases = gm_response(gm_response.collapse == 0,:);
    col_cases = gm_response(gm_response.collapse > 0,:);
    pt_diff = height(no_col_cases) -  height(col_cases);
    new_cases = col_cases(end,:);
    sa_orig = gm_response.sa_x(gm_response.scale == 1)
    new_sa = 0.25*sa_orig + col_cases.sa_x(end);
    new_cases.sa_x = new_sa;
    for s = 2:pt_diff
        new_sa = new_sa + 0.25*sa_orig;
        new_cases(s,:) = col_cases(end,:);
        new_cases.sa_x(s) = new_sa;
    end
    new_ida_table = [new_ida_table;no_col_cases;col_cases;new_cases];
end

%% For each component metric, quantify comparison logistic regression
% Collapse
attrs = {'sa_x', 'max_drift', 'drift_x', 'drift_z', 'gravity_load_lost_ratio','gravity_load_lost_ratio_alt','total_energy_ft_lbs','norm_energy_max','norm_energy_tot',...
        'num_adjacent_failure_any','num_adjacent_failure_any_frame','num_adjacent_failure_all', 'cols_walls_1_num_cp', 'cols_walls_1_percent_cp',...
        'cols_walls_1_max_cp','cols_walls_1_min_cp','cols_walls_1_mean_cp','cols_walls_1_range_cp','cols_walls_1_std_cp','cols_walls_1_cov_cp',... 
        'cols_walls_1_num_b_e','cols_walls_1_percent_b_e','cols_walls_1_max_b_e','cols_walls_1_mean_b_e'};
for a = 1:length(attrs)
    plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'ML Plots' '/' 'Collapse'];
    [LR.collapse.(attrs{a})] = fn_logistic_regression(max(new_ida_table.(attrs{a}),1e-6), min(new_ida_table.collapse,1), plot_dir, attrs{a});
end
    
    
% Save Logistic Regression Fits
save([read_dir filesep 'LR_data.mat'],'LR')

end

function [metrics] = fn_logistic_regression(max_component_response, collapse_target, plot_dir, x_label)
import plotting_tools.fn_format_and_save_plot

%% Logistic Regression
B = mnrfit(log(max_component_response),categorical(~collapse_target));
% B = glmfit(log(ida_table_full.sa_x),collapse_target,'binomial','link','probit');
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

max_val = max(max_component_response)*1.1;
x_points = linspace(min(max_component_response),max_val,1000);
z_predict = B(1) + log(x_points)*B(2);
prob_predict = 1 ./ (1+exp(-z_predict));

hold on
scatter(max_component_response,collapse_target,'DisplayName','Recorded Collapse')
plot(x_points,prob_predict,'DisplayName','Predicted Collapse')
xlabel(strrep(x_label,'_',' '))
ylabel('Probability of Collapse')
xlim([0,max_val])
fn_format_and_save_plot( plot_dir, x_label, 3, 1 )

% Change Metrics to table for outputs
metrics = struct2table(metrics);
end
