function [ adjusted_med_sa, beta, p_col_mce_adj, p_col_dbe_adj ] = fn_create_collapse_fragility(analysis, gm_set_table, ida_results, main_dir, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% import packages
import ida.*

%% Collect IDA data
id = 0;
id_missing = 0;
for gm = 1:height(gm_set_table)
    gm_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm))];
    sa_folders = dir([gm_dir filesep 'Sa_*']);
    for s = 1:length(sa_folders)
        % Load data
        outputs_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm)) '/' sa_folders(s).name];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
        story_file = [outputs_dir filesep 'story_analysis.mat'];
        if exist(outputs_file,'file') && exist(hinge_file,'file')
            load(outputs_file)
            load(hinge_file)
            load(story_file)
            if isfield(summary,'max_drift_x')
                id = id + 1;
                ida.id(id,1) = id;
                ida.eq_name{id,1} = gm_set_table.eq_name{gm};

                % X direction
                ida.sa_x(id,1) = summary.sa_x;
                ida.mce_ratio_x(id,1) = ida.sa_x(id,1)/ida_results.mce(1);
                ida.drift_x(id,1) = summary.max_drift_x;

                % z direction 
                if analysis.run_z_motion
                    ida.sa_z(id,1) = summary.sa_z;
                    ida.mce_ratio_z(id,1) = ida.sa_z(id,1)/ida_results.mce(2);
                    ida.drift_z(id,1) = summary.max_drift_z;
                end

                % max drift
                if analysis.run_z_motion
                    ida.max_drift(id,1) = max(summary.max_drift_x,summary.max_drift_z);
                else
                    ida.max_drift(id,1) = summary.max_drift_x;
                end

                % Collapse metrics
                ida.collapse(id,1) = summary.collapse;
                if summary.collapse > 0
                    ida.collapse_direction{id,1} = summary.collapse_direction;
                    ida.collapse_mech{id,1} = summary.collaspe_mech;
                    if strcmp(summary.collapse_direction,'x')
                        ida.collapse_x(id,1) = summary.collapse;
                        ida.collapse_z(id,1) = NaN;
                    elseif strcmp(summary.collapse_direction,'z')
                        ida.collapse_x(id,1) = NaN;
                        ida.collapse_z(id,1) = summary.collapse;
                    end
                    if contains(summary.collaspe_mech,'column')
                        ida.collapse_comps(id,1) = sum(hinge.b_ratio(strcmp(hinge.ele_type,'column')) >= 1);
                    elseif contains(summary.collaspe_mech,'wall')
                        ida.collapse_comps(id,1) = sum(hinge.e_ratio(strcmp(hinge.ele_type,'wall')) >= 1);
                    else
                        ida.collapse_comps(id,1) = sum(hinge.b_ratio >= 1);
                    end
                else
                    ida.collapse_direction{id,1} = 'NA';
                    ida.collapse_mech{id,1} = 'NA';
                    ida.collapse_x(id,1) = NaN;
                    ida.collapse_z(id,1) = NaN;
                    ida.collapse_comps(id,1) = NaN;
                end            
            else
                id_missing = id_missing + 1;
                missing_ida.eq_name{id_missing,1} = gm_set_table.eq_name{gm};
                missing_ida.sa(id_missing,1) = str2double(regexp(sa_folders(s).name,'(?<=_).+$','match'));
            end
        else
            id_missing = id_missing + 1;
            missing_ida.eq_name{id_missing,1} = gm_set_table.eq_name{gm};
            missing_ida.sa(id_missing,1) = str2double(regexp(sa_folders(s).name,'(?<=_).+$','match'));
        end
    end
end

ida_table = struct2table(ida);

% Remove all cases that failed to converge yet did not get far enough
failed_convergence = ida_table(ida_table.collapse == 5,:);
ida_table(ida_table.collapse == 5,:) = []; % filter out failed models

%% Define Collapse Fragility
num_gms = [];
num_collapse = [];
prob_col = [];
for sa = 1:length(analysis.sa_stripes)
    filt_sa = round(ida_table.sa_x,3) == round(analysis.sa_stripes(sa),3);
    num_gms(sa) = sum(filt_sa);
    num_collapse(sa) = sum(ida_table.collapse(filt_sa));
    prob_col(sa) = mean(ida_table.collapse(filt_sa));
end

% Fit lognormal distribution
[theta, beta] = fn_mle_pc(analysis.sa_stripes, num_gms, num_collapse);

% Find probability of Collapse at MCE and DBE
mu = log(theta);
p_col_mce = logncdf(ida_results.mce,mu,beta);
p_col_dbe = logncdf((2/3)*ida_results.mce,mu,beta);

% Adjust for P-695 values
ssf_table = readtable(['+ida' filesep 'p695_ssf_factor_sdc_d.csv']);
if analysis.run_z_motion
    factor_3D = 1.11;
else
    factor_3D = 1.0;
end
% mu_t_ew = 0.017 / 0.005; % rough estimate from pushover (may be slightly higher)
mu_t = 6; % assume slightly less than R for now, need to calc based on a pushover

% period_low = min(ssf_table.period - ida_results.period)
% period_high

ssf_table_filt = ssf_table{end,2:end};
mu_t_list = [1.0 1.1 1.5 2 3 4 6 8];
SSF = interp1(mu_t_list, ssf_table_filt, mu_t);
cmr = theta / ida_results.mce;
acmr = factor_3D * SSF * cmr;
adjusted_med_sa = acmr * ida_results.mce;

% Find probability of Collapse at MCE and DBE
mu = log(adjusted_med_sa);
beta_new = 0.6;
p_col_mce_adj = logncdf(ida_results.mce,mu,beta_new);
p_col_dbe_adj = logncdf((2/3)*ida_results.mce,mu,beta_new);

%% Save data

end