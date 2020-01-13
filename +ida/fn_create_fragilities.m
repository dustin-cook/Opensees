function [ ] = fn_create_fragilities(analysis, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
% params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
% mechs = { 'cols_1', 'walls_1', 'cols_walls_1'};
attrs = {'sa_x', 'max_drift', 'drift_x', 'gravity_load_lost_ratio',...
        'num_cp', 'percent_cp','max_cp','min_cp','mean_cp','range_cp','std_cp','cov_cp',... 
        'num_b','percent_b','min_b','max_b','mean_b'};
% params = {'b','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
params = {'b','io','ls','cp'};
frag_probs = [10 25 50 75 100];

%% Collect info for each ground motion (new fragilities)
ida_table = readtable([write_dir filesep 'ida_table.csv']);
gm_table = readtable([write_dir filesep 'gm_table.csv']);
gm_data.collapse = table;
gm_data.collapse_2 = table;
gm_data.gravity = table;
gm_data.collapse_x = table;
gm_data.collapse_z = table;

for gm = 1:height(gm_table)
    gm_response = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)),:);
    if ~isempty(gm_response)
        gm_response_no_col = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)) & ida_table.collapse == 0,:);
        gm_response_no_col.collapse(gm_response_no_col.sa_x == max(gm_response_no_col.sa_x) | gm_response_no_col.max_drift >= 0.06) = 1; % Flag just before collapse point
        jbc = gm_response_no_col;

        % Just Before Collapse
        gm_data.collapse.eq_name{gm} = gm_table.eq_name{gm};
        for a = 1:length(attrs)
            gm_data.collapse.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
        end
        
        % Post Process Collapse
        pp_collapse = jbc;
        pp_collapse.collapse(pp_collapse.max_drift >= 0.06) = 1;
        pp_collapse.collapse(pp_collapse.gravity_dcr > 0.4) = 1;
        gm_data.collapse_2.eq_name{gm} = gm_table.eq_name{gm};
        for a = 1:length(attrs)
            gm_data.collapse_2.(attrs{a})(gm) = min(pp_collapse.(attrs{a})(pp_collapse.collapse == 1));
        end

    %     % Full Gravity Exceedance
    %     attrs = {'sa_x', 'max_drift', 'drift_x', 'drift_z', 'gravity_load_lost_ratio','gravity_load_lost_ratio_alt','total_energy_ft_lbs','norm_energy_max','norm_energy_tot',...
    %             'num_adjacent_failure_any','num_adjacent_failure_any_frame','num_adjacent_failure_all', 'cols_walls_1_num_cp', 'cols_walls_1_percent_cp',...
    %             'cols_walls_1_max_cp','cols_walls_1_min_cp','cols_walls_1_mean_cp','cols_walls_1_range_cp','cols_walls_1_std_cp','cols_walls_1_cov_cp',... 
    %             'cols_walls_1_num_b_e','cols_walls_1_percent_b_e','cols_walls_1_min_b_e','cols_walls_1_max_b_e','cols_walls_1_mean_b_e'};
    %     gm_data.full_gravity.eq_name{gm} = gm_set_table.eq_name{gm};
    %     for a = 1:length(attrs)
    %         gm_data.full_gravity.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.gravity_dcr >= 1));
    %     end

    %     % Collapse x
    %     attrs = {'sa_x', 'max_drift', 'drift_x', 'drift_z', 'gravity_load_lost_ratio','gravity_load_lost_ratio_alt','total_energy_ft_lbs','norm_energy_ew',...
    %             'num_adjacent_failure_any','num_adjacent_failure_any_frame','num_adjacent_failure_all', 'cols_1_num_cp', 'cols_1_percent_cp',...
    %             'cols_1_max_cp','cols_1_min_cp','cols_1_mean_cp','cols_1_num_b','cols_1_percent_b','cols_1_max_b','cols_1_mean_b'};
    %     gm_data.collapse_x.eq_name{gm} = gm_set_table.eq_name{gm};
    %     for a = 1:length(attrs)
    %         if strcmp(collapse_dir,'x')
    %             gm_data.collapse_x.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
    %         else
    %             gm_data.collapse_x.(attrs{a})(gm) = NaN;
    %         end
    %     end

    %     % Collapse z
    %     attrs = {'sa_x', 'max_drift', 'drift_x', 'drift_z', 'gravity_load_lost_ratio','gravity_load_lost_ratio_alt','total_energy_ft_lbs','norm_energy_ns',...
    %             'num_adjacent_failure_any','num_adjacent_failure_any_frame','num_adjacent_failure_all', 'walls_1_num_cp', 'walls_1_percent_cp',...
    %             'walls_1_max_cp','walls_1_min_cp','walls_1_mean_cp','walls_1_num_e','walls_1_percent_e','walls_1_max_e','walls_1_mean_e'};
    %     gm_data.collapse_z.eq_name{gm} = gm_set_table.eq_name{gm};
    %     for a = 1:length(attrs)
    %         if strcmp(collapse_dir,'z')
    %             gm_data.collapse_z.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
    %         else
    %             gm_data.collapse_z.(attrs{a})(gm) = NaN;
    %         end
    %     end
    end
end

% Save Tabular Results as CSVs
save([write_dir filesep 'gm_data.mat'],'gm_data')

%% Create Fragility Curves based on Baker MLE (New Fragilities)
% collape and unnacceptable response
field_names = fieldnames(gm_data);
for i = 1:length(field_names)
    fld = field_names{i};
    for j = 2:width(gm_data.(fld))
        var_name = gm_data.(fld).Properties.VariableNames{j};
        [new_frag_curves.(fld).(var_name)] = fn_fit_fragility_MOM(gm_data.(fld).(var_name));
    end
end

% Save Frag Curve Data
save([write_dir filesep 'new_frag_curves.mat'],'new_frag_curves')

%% Collect info for each ground motion (traditional fragilities)
for gm = 1:height(gm_table)
    gm_response = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)),:);
    gm_table.sa_collapse(gm) = max([min(gm_response.sa_x(gm_response.collapse > 0)),NaN]);
    gm_table.sa_collapse_drift(gm) = max([min(gm_response.sa_x(gm_response.collapse == 1)),NaN]);
    gm_table.sa_collapse_drift_6(gm) = max([min(gm_response.sa_x(gm_response.max_drift >= 0.06)),NaN]);
    gm_table.sa_collapse_grav(gm) = max([min(gm_response.sa_x(gm_response.gravity_dcr > 0.4)),NaN]);
    gm_table.sa_collapse_2(gm) = max([min(gm_response.sa_x(gm_response.collapse > 0 | gm_response.max_drift >= 0.06 | gm_response.gravity_dcr > 0.4)),NaN]);
    gm_table.sa_collapse_convergence(gm) = max([min(gm_response.sa_x(gm_response.collapse == 3)),NaN]);
    gm_table.sa_UR_accept_15(gm) = max([min(gm_response.sa_x(gm_response.num_b_15 > 0)),NaN]);
    gm_table.sa_UR(gm) = max([min(gm_response.sa_x(gm_response.collapse == 1 | gm_response.collapse == 3 | gm_response.num_b_15 > 0)),NaN]);
    
    % 1% to 10% drift fragilities
    for d = 1:10 
%         gm_set_table.(['sa_drift_' num2str(d)])(gm) = max([min(gm_response.sa_x(gm_response.drift_x >= d/100 | gm_response.drift_z >= d/100)),NaN]);
        gm_table.(['sa_drift_' num2str(d)])(gm) = max([min(gm_response.sa_x(gm_response.drift_x >= d/100)),NaN]);
    end
    
    % grav load lost fragilities
    for f = 1:length(frag_probs)
        gm_table.(['sa_gravity_percent_lost_' num2str(frag_probs(f))])(gm) = max([min(gm_response.sa_x(gm_response.gravity_load_lost_ratio >= frag_probs(f)/100)),NaN]);
    end
    
%     % adjacent components
%     gm_set_table.sa_adjacent_component_any(gm) = max([min(gm_response.sa_x(gm_response.adjacent_failure_any == 1)),NaN]);
%     gm_set_table.sa_adjacent_component_any_frame(gm) = max([min(gm_response.sa_x(gm_response.adjacent_failure_any_frame == 1)),NaN]);
%     gm_set_table.sa_adjacent_component_all(gm) = max([min(gm_response.sa_x(gm_response.adjacent_failure_all == 1)),NaN]);
%     
%     % Energy
%     gm_table.collapse_energy(gm) = max([min(gm_response.total_energy_ft_lbs(gm_response.collapse == 1)),NaN]);
%     for f = 1:length(frag_probs)
%         gm_table.(['sa_collapse_energy_percent_' num2str(frag_probs(f))])(gm) = max([min(gm_response.sa_x(gm_response.total_energy_ft_lbs/gm_table.collapse_energy(gm) >= frag_probs(f)/100)),NaN]);
%     end
    
    % non directional component fragilities 
    for p = 1:length(params)
        gm_table.(['sa_first_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['num_' params{p}]) > 0)),NaN]);
        for pr = 1:length(frag_probs)
            gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
        end
    end
    
%     % X direction Curves
%     gm_set_table.sa_collapse_x(gm) = max([min(gm_response.sa_x(strcmp(gm_response.collapse_direction,'x'))),NaN]);
%     for p = 1:length(params)
%         if any(strcmp(gm_response.collapse_direction,'x'))
%             gm_set_table.(['sa_cols_1_first_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['cols_1_num_' params{p}]) > 0)),NaN]);
%         else
%             gm_set_table.(['sa_cols_1_first_' params{p}])(gm) = NaN;
%         end
%         for pr = 1:length(frag_probs)
%             if any(strcmp(gm_response.collapse_direction,'x'))
%                 gm_set_table.(['sa_cols_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['cols_1_percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
%             else
%                 gm_set_table.(['sa_cols_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = NaN;
%             end
%         end
%     end
%     
% 
%     % Z direction curves
%     if analysis.run_z_motion
%         gm_set_table.sa_collapse_z(gm) = max([min(gm_response.sa_x(strcmp(gm_response.collapse_direction,'z'))),NaN]);
%         for p = 1:length(params)
%             gm_set_table.(['sa_walls_1_first_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['walls_1_num_' params{p}]) > 0)),NaN]);
%             for pr = 1:length(frag_probs)
%                 gm_set_table.(['sa_walls_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = max([min(gm_response.sa_x(gm_response.(['walls_1_percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
%             end
%         end
%     end
end

% Save Tabular Results as CSVs
writetable(gm_table,[write_dir filesep 'gm_table.csv'])

%% Create Fragility Curves based on Baker MLE (traditional Fragilities)
% collape and unnacceptable response
[frag_curves.collapse] = fn_fit_fragility_MOM(gm_table.sa_collapse);
[frag_curves.collapse_drift] = fn_fit_fragility_MOM(gm_table.sa_collapse_drift);
[frag_curves.collapse_drift_6] = fn_fit_fragility_MOM(gm_table.sa_collapse_drift_6);
[frag_curves.collapse_grav] = fn_fit_fragility_MOM(gm_table.sa_collapse_grav);
[frag_curves.collapse_2] = fn_fit_fragility_MOM(gm_table.sa_collapse_2);
[frag_curves.collapse_convergence] = fn_fit_fragility_MOM(gm_table.sa_collapse_convergence);
[frag_curves.UR_accept_15] = fn_fit_fragility_MOM(gm_table.sa_UR_accept_15);
[frag_curves.UR] = fn_fit_fragility_MOM(gm_table.sa_UR);

% non directional component fragilities 
for p = 1:length(params)
    [frag_curves.(params{p})] = fn_multi_frag_curves(gm_table, params{p}, frag_probs, ida_table.num_comps(1));
end

% 1% to 5% drift fragilities
for d = 1:10 
    [frag_curves.drift.(['idr_' num2str(d)])] = fn_fit_fragility_MOM(gm_table.(['sa_drift_' num2str(d)]));
end

% grav load lost fragilities
for f = 1:length(frag_probs) 
    [frag_curves.gravity.(['percent_lost_' num2str(frag_probs(f))])] = fn_fit_fragility_MOM(gm_table.(['sa_gravity_percent_lost_' num2str(frag_probs(f))]));
end

% % adjacent components
% [frag_curves.adjacent_comp.any] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_any);
% [frag_curves.adjacent_comp.any_frame] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_any_frame);
% [frag_curves.adjacent_comp.all] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_all);

% % Collapse Energy
% for f = 1:length(frag_probs) 
%     [frag_curves.energy.(['percent_collapse_' num2str(frag_probs(f))])] = fn_fit_fragility_MOM(gm_table.(['sa_collapse_energy_percent_' num2str(frag_probs(f))]));
% end

% % X direction Curves
% [frag_curves.ew_collapse] = fn_fit_fragility_MOM(gm_set_table.sa_collapse_x);
% for p = 1:length(params)
%     [frag_curves.cols_1.(params{p})] = fn_multi_frag_curves(gm_set_table, 'cols_1', params{p}, frag_probs, ida_table.num_comps_cols_1(1));
% end
% 
% % Z direction curves
% if analysis.run_z_motion
%     [frag_curves.ns_collapse] = fn_fit_fragility_MOM(gm_set_table.sa_collapse_z);
%     for p = 1:length(params)
%         [frag_curves.walls_1.(params{p})] = fn_multi_frag_curves(gm_set_table, 'walls_1', params{p}, frag_probs, ida_table.num_comps_walls_1(1));
%     end
% end

% Save Frag Curve Data
save([write_dir filesep 'frag_curves.mat'],'frag_curves')

% Write outputs summary file
outputs.med_sa_collapse = frag_curves.collapse.theta;
outputs.med_sa_cp = frag_curves.cp.theta(1);
outputs.collapse_margin = outputs.med_sa_collapse/outputs.med_sa_cp;
outputs.num_comps = median(gm_data.collapse.num_cp);
outputs.percent_comps = median(gm_data.collapse.percent_cp);
outputs.drift = new_frag_curves.collapse.drift_x.theta;
writetable(struct2table(outputs),[write_dir filesep 'summary_outputs.csv'])
end

function [params] = fn_fit_fragility_MOM(limit_state_dist)
limit_state_dist(isnan(limit_state_dist)) = [];
limit_state_dist(limit_state_dist <= 0) = 1e-6;
if ~isempty(limit_state_dist)
    [pHat, ~] = lognfit(limit_state_dist);
    params.theta = exp(pHat(1));
    params.beta = pHat(2);
else
    params.theta = NaN;
    params.beta = NaN;
end
end

function [frag_curves] = fn_multi_frag_curves(gm_set_table, param, frag_probs, num_comp_mech)
frag_curves = table;
frag_curves.num_comp(1) = 1;
frag_curves.prct_mech(1) = round(1/num_comp_mech,3);
[fits] = fn_fit_fragility_MOM(gm_set_table.(['sa_' 'first_' param]));
frag_curves.theta(1) = fits.theta;
frag_curves.beta(1) = fits.beta;
for pr = 1:length(frag_probs)
    frag_curves.num_comp(pr+1) = ceil(num_comp_mech*frag_probs(pr)/100);
    frag_curves.prct_mech(pr+1) = frag_probs(pr)/100;
    [fits] = fn_fit_fragility_MOM(gm_set_table.(['sa_' num2str(frag_probs(pr)) '_percent_' param]));
    frag_curves.theta(pr+1) = fits.theta;
    frag_curves.beta(pr+1) = fits.beta;
end
end