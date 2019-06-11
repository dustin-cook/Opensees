%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
analysis.model_id = 11;
analysis.proceedure = 'NDP';
analysis.id = 38;
analysis.summit = 0;
analysis.run_ida = 0;
analysis.post_process_ida = 0;
analysis.plot_ida = 1;
analysis.gm_set = 'FEMA_far_field';

% IDA Inputs
hazard.curve.rp = [43, 72, 224, 475, 975, 2475, 4975];%[43, 72, 224, 475, 975, 2475, 4975];
hazard.curve.pga = [0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];%[0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];
analysis.collapse_drift = 0.06;

% Secondary options
analysis.dead_load = 1;
analysis.live_load = 1;
analysis.opensees_SP = 1;
analysis.type = 1;
analysis.nonlinear = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.03;
analysis.hinge_stiff_mod = 10;
analysis.run_eigen = 0;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 1;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

%% P695 Factors
% Spectra
Sms = 1.5;
Sm1 = 0.9;
Ts = 0.6;
p695_spectra.periods = [0.01:0.01:1, 1.1:0.1:3];
p695_spectra.Smt = Sms*ones(1,length(p695_spectra.periods));
p695_spectra.Smt(p695_spectra.periods > Ts) = Sm1 ./ p695_spectra.periods(p695_spectra.periods > Ts);

% Period Based Ductility
mu_t_ew = 0.016 / 0.005; % rough estimate from pushover (may be slightly higher)
mu_t_ns = 0.004 / 0.0017; % rough estimate from pushover (may be slightly higher)

% MCE
building_period.ew = 1.14;
building_period.ns = 0.44;
MCE_ew = interp1(p695_spectra.periods, p695_spectra.Smt, building_period.ew);
MCE_ns = interp1(p695_spectra.periods, p695_spectra.Smt, building_period.ns);

% SSF (based on table 7-1b)
SSF_ew = interp1([3,4], [1.275, 1.33],mu_t_ew);
SSF_ns = interp1([1.0 1.1 1.5 2 3 4 6 8], [1.00 1.05 1.1 1.13 1.18 1.22 1.28 1.33],mu_t_ew);

% Dispersion
beta_rtr = 0.4;
beta_tot = 0.6;

%% Initial Setup
% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Import packages
import plotting_tools.fn_format_and_save_plot

%% Define read and write directories
model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'model_data'];
tcl_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'opensees_data'];
asce41_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'asce_41_data'];

% Load in Model Tables
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
load([asce41_dir filesep 'element_analysis.mat']);

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);
IDA_scale_factors = hazard.curve.pga ./ gm_median_pga;

%% Run Opensees Models
if analysis.run_ida || analysis.post_process_ida
parpool; % Set up Parallel Workers
for i = 1:length(IDA_scale_factors)
    error_count = 0;
    scale_factor = IDA_scale_factors(i);
    parfor gms = 1:height(gm_set_table)
        % Run Opensees
        if analysis.run_ida
            % Suppress MATLAB warnings
            warning('off','all')
            fprintf('Running Scale Factor %4.2f for Ground Motion ID: %i-%i \n\n', scale_factor, gm_set_table.set_id(gms), gm_set_table.pair(gms))
            [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, building_period, tcl_dir);
            if exit_status == 1
                error_count = error_count + 1;
            end
        else
            exit_status = 0;
        end

        if analysis.post_process_ida && exit_status ~= 1
            fprintf('Postprocessing Opensees Ouputs\n')
            fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor)
        end
        fprintf('\n')
    end
    fprintf('%i Failed GMs for Scale Factor %4.2f \n\n', error_count, scale_factor)
end
delete(gcp('nocreate')) % End Parallel Process
end



%% Plot Results
if analysis.plot_ida
    col_ids = [1 3 5 7 9 11 12 14 16 18 20 22 23 25 27 29 31 33 34 36 38 40 42 44];
    % Plot Results
    plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'IDA Plots'];
    id = 0;
    id_missing = 0;
    for i = 1:length(IDA_scale_factors)
        for gms = 1:height(gm_set_table)   
            % Load data
            outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(IDA_scale_factors(i)) '/' 'GM_' num2str(gm_set_table.set_id(gms)) '_' num2str(gm_set_table.pair(gms))];
            outputs_file = [outputs_dir filesep 'summary_results.mat'];
            if exist(outputs_file,'file')
                id = id + 1;
                load(outputs_file)
                load([outputs_dir filesep 'hinge_analysis.mat'])
                ida.id(id,1) = id;
                ida.scale(id,1) = IDA_scale_factors(i);
                ida.gm_name{id,1} = gm_set_table.eq_name{gms};
                ida.sa_x(id,1) = summary.sa_x;
                ida.sa_z(id,1) = summary.sa_z;
                ida.drift_x(id,1) = summary.max_drift_x;
                ida.drift_z(id,1) = summary.max_drift_z;
                ida.collapse(id,1) = summary.collapse;
                col_base_filter = hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'primary') & ismember(hinge.element_id,col_ids);
                col_b_ratios = hinge.b_ratio(col_base_filter);
                ida.num_b(id,1) = sum(col_b_ratios >= 1.0);
                ida.percent_b(id,1) = ida.num_b(id,1) / length(col_b_ratios);
                ida.num_15b(id,1) = sum(col_b_ratios >= 1.5);
                ida.percent_15b(id,1) = ida.num_15b(id,1) / length(col_b_ratios);
                col_io_ratios = hinge.io_ratio(col_base_filter);
                ida.num_io(id,1) = sum(col_io_ratios >= 1.0);
                ida.percent_io(id,1) = ida.num_io(id,1) / length(col_io_ratios);
                col_ls_ratios = hinge.ls_ratio(col_base_filter);
                ida.num_ls(id,1) = sum(col_ls_ratios >= 1.0);
                ida.percent_ls(id,1) = ida.num_ls(id,1) / length(col_ls_ratios);                
                col_cp_ratios = hinge.cp_ratio(col_base_filter);
                ida.num_cp(id,1) = sum(col_cp_ratios >= 1.0);
                ida.percent_cp(id,1) = ida.num_cp(id,1) / length(col_cp_ratios);
                col_euro_ratios = hinge.max_deform(col_base_filter) ./ hinge.euro_th_NC(col_base_filter); % NEED TO UPDATE THIS FOR NEXT RUN
                ida.num_euro(id,1) = sum(col_euro_ratios >= 1.0);
                ida.percent_euro(id,1) = ida.num_euro(id,1) / length(col_euro_ratios);
            else
                id_missing = id_missing + 1;
                missing_ida.scale(id_missing,1) = IDA_scale_factors(i);
                missing_ida.gm_set_id(id_missing,1) = gm_set_table.set_id(gms);
                missing_ida.gm_set_pair_id(id_missing,1) = gm_set_table.pair(gms);
            end
        end
    end

    %% Plot IDA curves
    fprintf('Saving IDA Summary Data and Figures to Directory: %s',plot_dir)
    hold on
    for gms = 1:height(gm_set_table)
        plot(ida.drift_x(strcmp(ida.gm_name,gm_set_table.eq_name{gms})),ida.sa_x(strcmp(ida.gm_name,gm_set_table.eq_name{gms})))
    end
    xlabel('Max Drift')
    ylabel('Sa(T_1) (g)')
    fn_format_and_save_plot( plot_dir, 'IDA Plot EW Frame Direction', 2 )
    hold on
    for gms = 1:height(gm_set_table)
        plot(ida.drift_z(strcmp(ida.gm_name,gm_set_table.eq_name{gms})),ida.sa_z(strcmp(ida.gm_name,gm_set_table.eq_name{gms})))
    end
    xlabel('Max Drift')
    ylabel('Sa(T_1) (g)')
    fn_format_and_save_plot( plot_dir, 'IDA Plot NS Wall Direction', 2 )

    %% Save Tabular Results as CSV
    ida_table = struct2table(ida);
    writetable(ida_table,[plot_dir filesep 'ida_results.csv'])
    missing_ida_table = struct2table(missing_ida);
    writetable(missing_ida_table,[plot_dir filesep 'missing_ida_results.csv'])
    
    %% Collect frag curve info
    % Remove all cases that failed to converge yet did not get far enough
    ida_table(ida_table.collapse == 5,:) = [];
    
    % Collect Stripe info
    for i = 1:length(IDA_scale_factors)
        % Regular Collapse
        collapse_ids = ida_table.collapse(ida_table.scale == IDA_scale_factors(i));
        frag.num_gms(i,1) = length(collapse_ids);
        frag.num_collapse(i,1) = sum(collapse_ids ~= 0);
        frag.prob_collapse(i,1) = frag.num_collapse(i,1) / frag.num_gms(i,1);
        
        % B collapse
        percent_b = ida_table.percent_b(ida_table.scale == IDA_scale_factors(i));
        num_b = ida_table.num_b(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_b(i,1) = sum(num_b >= 1) / length(num_b);
        frag.prob_collapse_b_10(i,1) = sum(percent_b >= 0.1) / length(percent_b);
        frag.prob_collapse_b_25(i,1) = sum(percent_b >= 0.25) / length(percent_b);
        frag.prob_collapse_b_50(i,1) = sum(percent_b >= 0.5) / length(percent_b);
        frag.prob_collapse_b_100(i,1) = sum(percent_b >= 1.0) / length(percent_b);
        
        % 1.5B collapse
        percent_15b = ida_table.percent_15b(ida_table.scale == IDA_scale_factors(i));
        num_15b = ida_table.num_15b(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_15b(i,1) = sum(num_15b >= 1) / length(num_15b);
        frag.prob_collapse_15b_10(i,1) = sum(percent_15b >= 0.1) / length(percent_15b);
        frag.prob_collapse_15b_25(i,1) = sum(percent_15b >= 0.25) / length(percent_15b);
        frag.prob_collapse_15b_50(i,1) = sum(percent_15b >= 0.5) / length(percent_15b);
        frag.prob_collapse_15b_100(i,1) = sum(percent_15b >= 1.0) / length(percent_15b);
        
        % Acceptance Criteria
        num_io = ida_table.num_io(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_io(i,1) = sum(num_io >= 1) / length(num_io);
        num_ls = ida_table.num_ls(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_ls(i,1) = sum(num_ls >= 1) / length(num_ls);
        num_cp = ida_table.num_cp(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_cp(i,1) = sum(num_cp >= 1) / length(num_cp);
        num_euro = ida_table.num_euro(ida_table.scale == IDA_scale_factors(i));
        frag.prob_collapse_first_euro(i,1) = sum(num_euro >= 1) / length(num_euro);
        
        % Spectral Accelration
        frag.med_sa_ew(i,1) = median(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)));
        frag.med_sa_ns(i,1) = median(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)));
    end
    
    %% Calculate Post Fragulity Curve P695 factors
    % Collapse Median Sa
    Sa_med_col_ew = max(frag.med_sa_ew(frag.prob_collapse <= 0.5));
    Sa_med_col_ns = max(frag.med_sa_ns(frag.prob_collapse <= 0.5));
    
    % Collapse Margin Ratio
    CMR_ew = Sa_med_col_ew / MCE_ew;
    CMR_ns = Sa_med_col_ns / MCE_ns;
    
    % Adjust for SSF and 3D
    ACMR_ew = 1.2*SSF_ew*CMR_ew;
    ACMR_ns = 1.2*SSF_ns*CMR_ns;
    
    % Create Fragility curves based on P695
    x_points = 0:0.01:2;
    p695_frag_curve_ew = logncdf(x_points,log(Sa_med_col_ew),beta_rtr);
    p695_frag_curve_ns = logncdf(x_points,log(Sa_med_col_ns),beta_rtr);
    
    % Create Fragility Curves based on Baker MLE
    [ML_theta_ew, MLE_beta_ew] = fn_mle_pc(frag.med_sa_ew, frag.num_gms, frag.num_collapse);
    MLE_frag_curve_ew = logncdf(x_points,log(ML_theta_ew),MLE_beta_ew);
    [ML_theta_ns, MLE_beta_ns] = fn_mle_pc(frag.med_sa_ns, frag.num_gms, frag.num_collapse);
    MLE_frag_curve_ns = logncdf(x_points,log(ML_theta_ns),MLE_beta_ns);
    
    %% Plot Frag Curves
    figure
    hold on
    title('P695 Collapse Fragility: 6% Drift') 
    scatter(frag.med_sa_ew,frag.prob_collapse,'b','filled','DisplayName','EW Frame')
    plot(x_points,p695_frag_curve_ew,'b','DisplayName','EW Frame - P695 Curve')
    plot(x_points,MLE_frag_curve_ew,'--b','DisplayName','EW Frame - MLE Curve')
    scatter(frag.med_sa_ns,frag.prob_collapse,'k','filled','DisplayName','NS Wall')
    plot(x_points,p695_frag_curve_ns,'k','DisplayName','NS Wall - P695 Curve')
    plot(x_points,MLE_frag_curve_ns,'--k','DisplayName','NS Wall - MLE Curve')
    xlabel('Median Sa (g)')
    ylabel('P[Collapse]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility', 5 )
    
    figure
    hold on
    title('Frist Story Column Base B-Value Fragility: EW') 
    plot(frag.med_sa_ew,frag.prob_collapse_first_b,'DisplayName','First Column Base')
    plot(frag.med_sa_ew,frag.prob_collapse_b_10,'DisplayName','10% of Column Bases')
    plot(frag.med_sa_ew,frag.prob_collapse_b_25,'DisplayName','25% of Columns')
    plot(frag.med_sa_ew,frag.prob_collapse_b_50,'DisplayName','50% of Columns')
    plot(frag.med_sa_ew,frag.prob_collapse_b_100,'DisplayName','100% of Columns')
    xlabel('Median Sa (g)')
    ylabel('P[Exceeding Column Base B-Value]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility - b ratio - EW', 5 )
    
    figure
    hold on
    title('Frist Story Column Base 1.5xB-Value Fragility: EW')
    plot(frag.med_sa_ew,frag.prob_collapse_first_15b,'DisplayName','First Column Base')
    plot(frag.med_sa_ew,frag.prob_collapse_15b_10,'DisplayName','10% of Column Bases')
    plot(frag.med_sa_ew,frag.prob_collapse_15b_25,'DisplayName','25% of Column Bases')
    plot(frag.med_sa_ew,frag.prob_collapse_15b_50,'DisplayName','50% of Column Bases')
    plot(frag.med_sa_ew,frag.prob_collapse_15b_100,'DisplayName','100% of Column Bases')
    xlabel('Median Sa (g)')
    ylabel('P[Exceeding Column Base 1.5xB-Value]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility - 1_5b ratio - EW', 5 )
    
    figure
    hold on
    title('Frist Story Column Base Acceptance Criteria: EW') 
    plot(frag.med_sa_ew,frag.prob_collapse_first_io,'DisplayName','First IO Column')
    plot(frag.med_sa_ew,frag.prob_collapse_first_ls,'DisplayName','First LS Column')
    plot(frag.med_sa_ew,frag.prob_collapse_first_cp,'DisplayName','First CP Column')
    plot(frag.med_sa_ew,frag.prob_collapse_first_euro,'DisplayName','First Euro Column')
    xlabel('Median Sa (g)')
    ylabel('P[Exceeding Column Base Accepatance Criteria]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility - Acceptance Criteria - EW', 5 )
    
    % Save Collapse Fragilty Table
    

end
