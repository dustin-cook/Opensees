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
analysis.id = 3;
analysis.summit = 0;
analysis.run_ida = 1;
analysis.post_process_ida = 0;
analysis.plot_ida = 0;
analysis.gm_set = 'FEMA_far_field';

% IDA Inputs
% hazard.curve.rp = [22, 35, 64, 108, 144];
% hazard.curve.pga = [0.128, 0.192, 0.288, 0.376, 0.425];
hazard.curve.rp = [22, 35, 43, 64, 72, 108, 144, 224, 475, 975, 2475, 4975];
hazard.curve.pga = [0.128, 0.192, 0.224, 0.288, 0.308, 0.376, 0.425 0.502, 0.635, 0.766, 0.946, 1.082];
% hazard.curve.rp = [43, 72, 224, 475, 975, 2475, 4975];
% hazard.curve.pga = [0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];
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
% p695_spectra.periods = [0.01:0.01:1, 1.1:0.1:3];
% p695_spectra.Smt = Sms*ones(1,length(p695_spectra.periods));
% p695_spectra.Smt(p695_spectra.periods > Ts) = Sm1 ./ p695_spectra.periods(p695_spectra.periods > Ts);

% Period Based Ductility
mu_t_ew = 0.016 / 0.005; % rough estimate from pushover (may be slightly higher)
mu_t_ns = 0.004 / 0.0017; % rough estimate from pushover (may be slightly higher)

% MCE
ida_results.direction = {'EW'; 'NS'};
ida_results.period = [1.14; 0.44];
ida_results.spectra = [0.35; 0.77];
ida_results.mce = [0.53; 1.36]; % MCE Max from SP3, fixed to this site and model period

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
% parpool; % Set up Parallel Workers
for i = 1:length(IDA_scale_factors)
    error_count = 0;
    scale_factor = IDA_scale_factors(i);
    analysis.ground_motion_scale_factor = scale_factor;
    run_ida = analysis.run_ida;
    for gms = 1:height(gm_set_table)
        % Run Opensees
        if run_ida
            % Suppress MATLAB warnings
            warning('off','all')
            fprintf('Running Scale Factor %4.2f for Ground Motion ID: %i-%i \n\n', scale_factor, gm_set_table.set_id(gms), gm_set_table.pair(gms))
            [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, ida_results.period, tcl_dir);
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
% delete(gcp('nocreate')) % End Parallel Process
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
            hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
            if exist(outputs_file,'file') && exist(hinge_file,'file')
                id = id + 1;
                load(outputs_file)
                load(hinge_file)
                ida.id(id,1) = id;
                ida.scale(id,1) = IDA_scale_factors(i);
                ida.gm_name{id,1} = gm_set_table.eq_name{gms};
                ida.sa_x(id,1) = summary.sa_x;
                ida.sa_z(id,1) = summary.sa_z;
                ida.drift_x(id,1) = summary.max_drift_x;
                ida.drift_z(id,1) = summary.max_drift_z;
                if ida.drift_x(id,1) > 0.1 || ida.drift_z(id,1) > 0.1
                    ida.collapse(id,1) = 5; % Omit 1 case of extreme drift that doesnt make sense
                else
                    ida.collapse(id,1) = summary.collapse;
                end
                if isfield(summary,'collapse_direction')
                    ida.collapse_direction{id,1} = summary.collapse_direction;
                end
                col_base_filter = hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'primary') & ismember(hinge.element_id,col_ids);
                col_hinges = hinge(col_base_filter,:);
                col_b_ratios = col_hinges.b_ratio;
                ida.b_cols{id,1} = col_hinges.node_1(col_b_ratios >= 1.0);
                if summary.collapse == 3 || summary.collapse == 1
                    ida.num_b(id,1) = length(col_b_ratios);
                else
                    ida.num_b(id,1) = sum(col_b_ratios >= 1.0);
                end
                ida.percent_b(id,1) = ida.num_b(id,1) / length(col_b_ratios);
                col_io_ratios = hinge.io_ratio(col_base_filter);
                ida.io_cols{id,1} = col_hinges.node_1(col_io_ratios >= 1.0);
                if summary.collapse == 3 || summary.collapse == 1
                    ida.num_io(id,1) = length(col_io_ratios);
                else
                    ida.num_io(id,1) = sum(col_io_ratios >= 1.0);
                end
                ida.percent_io(id,1) = ida.num_io(id,1) / length(col_io_ratios);
                col_ls_ratios = hinge.ls_ratio(col_base_filter);
                ida.ls_cols{id,1} = col_hinges.node_1(col_ls_ratios >= 1.0);
                if summary.collapse == 3 || summary.collapse == 1
                    ida.num_ls(id,1) = length(col_ls_ratios);
                else
                    ida.num_ls(id,1) = sum(col_ls_ratios >= 1.0);
                end
                ida.percent_ls(id,1) = ida.num_ls(id,1) / length(col_ls_ratios);                
                col_cp_ratios = hinge.cp_ratio(col_base_filter);
                ida.cp_cols{id,1} = col_hinges.node_1(col_cp_ratios >= 1.0);
                if summary.collapse == 3 || summary.collapse == 1
                    ida.num_cp(id,1) = length(col_cp_ratios);
                else
                    ida.num_cp(id,1) = sum(col_cp_ratios >= 1.0);
                end
                ida.percent_cp(id,1) = ida.num_cp(id,1) / length(col_cp_ratios);
                col_euro_ratios = hinge.euro_V_NC_ratio;
                ida.euro_cols{id,1} = col_hinges.node_1(col_euro_ratios >= 1.0);
                if summary.collapse == 3 || summary.collapse == 1
                    ida.num_euro(id,1) = length(col_euro_ratios);
                else
                    ida.num_euro(id,1) = sum(col_euro_ratios >= 1.0);
                end
                ida.percent_euro(id,1) = ida.num_euro(id,1) / length(col_euro_ratios);
            else
                id_missing = id_missing + 1;
                missing_ida.scale(id_missing,1) = IDA_scale_factors(i);
                missing_ida.gm_set_id(id_missing,1) = gm_set_table.set_id(gms);
                missing_ida.gm_set_pair_id(id_missing,1) = gm_set_table.pair(gms);
            end
        end
    end

    % filter non_collapse 
    % Remove all cases that failed to converge yet did not get far enough
    ida_table = struct2table(ida);
    failed_convergence = ida_table(ida_table.collapse == 5,:);
    ida_table(ida_table.collapse == 5,:) = [];
    
    % Save Tabular Results as CSVs
    writetable(ida_table,[plot_dir filesep 'ida_table.csv'])
    writetable(struct2table(missing_ida),[plot_dir filesep 'idas_missing.csv'])
    writetable(failed_convergence,[plot_dir filesep 'idas_failed_convergence.csv'])
    
    %% Calculate Median Drifts and Accels
    for i = 1:length(IDA_scale_factors)
        % Spectral Accelration
        frag.med_sa_ew(i,1) = exp(mean(log(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)))));
        frag.p_15_sa_ew(i,1) = prctile(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)),15);
        frag.p_85_sa_ew(i,1) = prctile(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)),85);
        frag.med_sa_ns(i,1) = exp(mean(log(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)))));
        frag.p_15_sa_ns(i,1) = prctile(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)),15);
        frag.p_85_sa_ns(i,1) = prctile(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)),85);

        % Drift
        frag.mean_idr_ew(i,1) = mean(ida_table.drift_x(ida_table.scale == IDA_scale_factors(i)));
        frag.mean_idr_ns(i,1) = exp(mean(log(ida_table.drift_z(ida_table.scale == IDA_scale_factors(i)))));
    end
    
    %% Plot IDA curves
    fprintf('Saving IDA Summary Data and Figures to Directory: %s',plot_dir)
    hold on
    for gms = 1:height(gm_set_table)
        ida_plt = plot(ida_table.drift_x(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),ida_table.sa_x(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),'color',[0.75 0.75 0.75]);
        set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    plot(frag.mean_idr_ew,frag.med_sa_ew,'b','lineWidth',1.5,'DisplayName','Mean Drift')
    plot(frag.mean_idr_ew,frag.p_15_sa_ew,'--b','lineWidth',1.5,'DisplayName','15th Percentile')
    plot(frag.mean_idr_ew,frag.p_85_sa_ew,'--b','lineWidth',1.5,'DisplayName','85th Percentile')
    plot([0,0.08],[ida_results.spectra(1),ida_results.spectra(1)],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
%     set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    xlim([0 0.08])
    xlabel('Max Drift')
    ylabel('Sa(T_1=1.14s) (g)')
    fn_format_and_save_plot( plot_dir, 'IDA Plot EW Frame Direction', 1 )
    hold on
    for gms = 1:height(gm_set_table)
        ida_plt = plot(ida_table.drift_z(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),ida_table.sa_z(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),'color',[0.75 0.75 0.75]);
        set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    plot(frag.mean_idr_ns,frag.med_sa_ns,'b','lineWidth',1.5,'DisplayName','Mean Drift')
    plot(frag.mean_idr_ns,frag.p_15_sa_ns,'--b','lineWidth',1.5,'DisplayName','15th Percentile')
    plot(frag.mean_idr_ns,frag.p_85_sa_ns,'--b','lineWidth',1.5,'DisplayName','85th Percentile')
    plot([0,0.08],[ida_results.spectra(2),ida_results.spectra(2)],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
%     set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    xlim([0 0.08])
    xlabel('Max Drift')
    ylabel('Sa(T_1=0.44s) (g)')
    fn_format_and_save_plot( plot_dir, 'IDA Plot NS Wall Direction', 1 )
    
    %% Collect frag curve info
    % Collect Stripe info
    frag_b = [];
    frag_io = [];
    frag_ls = [];
    frag_cp = [];
    frag_euro = [];
    frag_b_drift = [];
    for i = 1:length(IDA_scale_factors)
        % Regular Collapse
        collapse_ids = ida_table.collapse(ida_table.scale == IDA_scale_factors(i));
        frag.num_gms(i,1) = length(collapse_ids);
        frag.num_gms_drift(i,1) = sum(ida_table.scale == IDA_scale_factors(i) & ida_table.collapse ~= 3);
        frag.num_collapse_drift(i,1) = sum(collapse_ids == 1);
        frag.num_collapse_convergence(i,1) = sum(collapse_ids == 3);
        frag.num_collapse_tot(i,1) = frag.num_collapse_drift(i,1) + frag.num_collapse_convergence(i,1);
        frag.prob_collapse_drift(i,1) = frag.num_collapse_drift(i,1) / frag.num_gms(i,1);
        frag.prob_collapse_convergence(i,1) = frag.num_collapse_convergence(i,1) / frag.num_gms(i,1);
        frag.prob_collapse_tot(i,1) = frag.num_collapse_tot(i,1) / frag.num_gms(i,1);
        
        % B collapse
        percent_b = ida_table.percent_b(ida_table.scale == IDA_scale_factors(i));
        num_b = ida_table.num_b(ida_table.scale == IDA_scale_factors(i));
        [frag_b] = fn_exceedance_values(frag_b, i, num_b, percent_b);
        percent_b_drift = ida_table.percent_b(ida_table.scale == IDA_scale_factors(i) & ida_table.collapse ~= 3);
        num_b_drift = ida_table.num_b(ida_table.scale == IDA_scale_factors(i) & ida_table.collapse ~= 3);
        [frag_b_drift] = fn_exceedance_values(frag_b_drift, i, num_b_drift, percent_b_drift);
        
        % IO
        percent_io = ida_table.percent_io(ida_table.scale == IDA_scale_factors(i));
        num_io = ida_table.num_io(ida_table.scale == IDA_scale_factors(i));
        [frag_io] = fn_exceedance_values(frag_io, i, num_io, percent_io);
        
        % LS
        percent_ls = ida_table.percent_ls(ida_table.scale == IDA_scale_factors(i));
        num_ls = ida_table.num_ls(ida_table.scale == IDA_scale_factors(i));
        [frag_ls] = fn_exceedance_values(frag_ls, i, num_ls, percent_ls);
        
        % CP
        percent_cp = ida_table.percent_cp(ida_table.scale == IDA_scale_factors(i));
        num_cp = ida_table.num_cp(ida_table.scale == IDA_scale_factors(i));
        [frag_cp] = fn_exceedance_values(frag_cp, i, num_cp, percent_cp);
        
        % Eurocode
        percent_euro = ida_table.percent_euro(ida_table.scale == IDA_scale_factors(i));
        num_euro = ida_table.num_euro(ida_table.scale == IDA_scale_factors(i));
        [frag_euro] = fn_exceedance_values(frag_euro, i, num_euro, percent_euro);
    end
    
    %% Calculate Post Fragulity Curve P695 factors
    % Collapse Median Sa
    ida_results.sa_med_col(1,1) = max(frag.med_sa_ew(frag.prob_collapse_tot <= 0.5));
    ida_results.sa_med_col(2,1) = max(frag.med_sa_ns(frag.prob_collapse_tot <= 0.5));
    
    % Collapse Margin Ratio
    ida_results.cmr(1,1) = ida_results.sa_med_col(1) / ida_results.mce(1);
    ida_results.cmr(2,1) = ida_results.sa_med_col(2) / ida_results.mce(2);
    
    % Adjust for SSF and 3D
    ida_results.acmr(1,1) = 1.2*SSF_ew*ida_results.cmr(1);
    ida_results.acmr(2,1) = 1.2*SSF_ns*ida_results.cmr(2);
    median_adjustment(1) = ida_results.sa_med_col(1)*(1.2*SSF_ew - 1);
    median_adjustment(2) = ida_results.sa_med_col(2)*(1.2*SSF_ns - 1);
    
    % Create Fragility curves based on P695
    x_points = 0:0.01:2;
    x_points_adjused_x = x_points + median_adjustment(1);
    x_points_adjused_z = x_points + median_adjustment(2);
%     p695_frag_curve_ew = logncdf(x_points,log(Sa_med_col_ew),beta_rtr);
%     p695_frag_curve_ns = logncdf(x_points,log(Sa_med_col_ns),beta_rtr);
    
    % Create Fragility Curves based on Baker MLE
    [col_frag_curve_ew_drift, ~] = fn_calc_frag_curve(x_points, frag.med_sa_ew, frag.num_gms_drift, frag.num_collapse_drift);
    [col_frag_curve_ns_drift, ~] = fn_calc_frag_curve(x_points, frag.med_sa_ns, frag.num_gms_drift, frag.num_collapse_drift);
    [col_frag_curve_ew_convergence, ~] = fn_calc_frag_curve(x_points, frag.med_sa_ew, frag.num_gms, frag.num_collapse_convergence);
    [col_frag_curve_ns_convergence, ~] = fn_calc_frag_curve(x_points, frag.med_sa_ns, frag.num_gms, frag.num_collapse_convergence);
    [col_frag_curve_ew_tot, col_frag_curve_full_ew_tot] = fn_calc_frag_curve(x_points, frag.med_sa_ew, frag.num_gms, frag.num_collapse_tot);
    [col_frag_curve_ns_tot, col_frag_curve_full_ns_tot] = fn_calc_frag_curve(x_points, frag.med_sa_ns, frag.num_gms, frag.num_collapse_tot);
    [frag_b] = fn_multi_frag_curves(x_points, frag_b, frag.med_sa_ew, frag.num_gms);
    [frag_b_drift] = fn_multi_frag_curves(x_points, frag_b_drift, frag.med_sa_ew, frag.num_gms_drift);
    [frag_io] = fn_multi_frag_curves(x_points, frag_io, frag.med_sa_ew, frag.num_gms);
    [frag_ls] = fn_multi_frag_curves(x_points, frag_ls, frag.med_sa_ew, frag.num_gms);
    [frag_cp] = fn_multi_frag_curves(x_points, frag_cp, frag.med_sa_ew, frag.num_gms);
    [frag_euro] = fn_multi_frag_curves(x_points, frag_euro, frag.med_sa_ew, frag.num_gms);
    
    % Adjust Collapse curves based on  
    
    
    %% Save IDA results as table
    [~, idx] = min(abs(x_points_adjused_x - ida_results.mce(1)));
    ida_results.p_col_mce(1,1) = col_frag_curve_ew_tot(idx);
    [~, idx] = min(abs(x_points_adjused_z - ida_results.mce(2)));
    ida_results.p_col_mce(2,1) = col_frag_curve_ns_tot(idx);
    ida_results_table = struct2table(ida_results);
    writetable(ida_results_table,[plot_dir filesep 'ida_results.csv'])
    
    %% Plot Frag Curves    
    % Total Collapse with discrete scatter
    figure
    hold on
    title('P695 Collapse Fragility') 
    sct = scatter(frag.med_sa_ew,frag.prob_collapse_tot,'b','filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,col_frag_curve_ew_tot,'b','DisplayName','EW Frame')
    plot([ida_results.spectra(1),ida_results.spectra(1)],[0,1],':b','lineWidth',1.5,'DisplayName','EW Frame - ICSB Level')
    sct = scatter(frag.med_sa_ns,frag.prob_collapse_tot,'k','filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,col_frag_curve_ns_tot,'k','DisplayName','NS Wall')
    plot([ida_results.spectra(2),ida_results.spectra(2)],[0,1],':k','lineWidth',1.5,'DisplayName','NS Wall - ICSB Level')
    xlabel('Sa(T_1) (g)')
    ylabel('P[Collapse]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility', 5 )
    
    % Total Collapse with 695 adjustments
    figure
    hold on
    title('P695 Collapse Fragility') 
    plot(x_points,col_frag_curve_ew_tot,'b','DisplayName','EW Frame')
    plot(x_points_adjused_x,col_frag_curve_ew_tot,'--b','DisplayName','EW Frame - Adjusted')
    plot(x_points_adjused_x,col_frag_curve_full_ew_tot,':b','DisplayName','EW Frame - Full Uncertainty')
    plot(x_points,col_frag_curve_ns_tot,'k','DisplayName','NS Wall')
    plot(x_points_adjused_z,col_frag_curve_ns_tot,'--k','DisplayName','NS Wall - Adjusted')
    plot(x_points_adjused_z,col_frag_curve_full_ns_tot,':k','DisplayName','NS Wall - Full Uncertainty')
    xlabel('Sa(T_1) (g)')
    ylabel('P[Collapse]')
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility with 695 Uncertainty', 5 )
    
    % Total v 6% v Convergence
    figure
    hold on
    title('P695 Collapse Fragility: EW Frame') 
    plot(frag.med_sa_ew,frag.prob_collapse_tot,'k','DisplayName','Total');
    plot(frag.med_sa_ew,frag.prob_collapse_drift,'b','DisplayName','6% Drift');
    plot(frag.med_sa_ew,frag.prob_collapse_convergence,'r','DisplayName','Nonconvergence');
    xlabel('Sa(T_1=1.14) (g)')
    ylabel('P[Collapse]')
    xlim([0 1])
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility Breakdown - EW', 5 )
    
    % Column B Value
    figure
    hold on
    title('First Story Column Base B-Value Fragility: EW')
    plot_name = 'Collapse Fragility - b ratio - EW';
    fn_plot_frag_curves(x_points, frag_b, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)

    % Column B Value - just drift
    figure
    hold on
    title('First Story Column Base B-Value Fragility: EW')
    plot_name = 'Collapse Fragility - b ratio - EW - Ignoring Convergence';
    fn_plot_frag_curves(x_points, frag_b_drift, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)
    
    % Column IO
    figure
    hold on
    title('First Story Column Base IO Fragility: EW')
    plot_name = 'Collapse Fragility - io - EW';
    fn_plot_frag_curves(x_points, frag_io, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)
    
    % Column LS
    figure
    hold on
    title('First Story Column Base LS Fragility: EW')
    plot_name = 'Collapse Fragility - ls - EW';
    fn_plot_frag_curves(x_points, frag_ls, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)
    
    % Column CP
    figure
    hold on
    title('First Story Column Base CP Fragility: EW')
    plot_name = 'Collapse Fragility - cp - EW';
    fn_plot_frag_curves(x_points, frag_cp, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)
    
    % Column Eurocode
    figure
    hold on
    title('First Story Column Base Eurocode Fragility: EW')
    plot_name = 'Collapse Fragility - euro - EW';
    fn_plot_frag_curves(x_points, frag_euro, frag.med_sa_ew, col_frag_curve_ew_tot, col_frag_curve_ew_drift, ida_results.spectra(1), plot_name, plot_dir)
    
    % Combined acceptance plots
%     color4 = [166,97,26
%       223,194,125
%       128,205,193
%       1,133,113]/255;
  
      matlab_colors = [0 0.4470 0.7410;
              0.85 0.325 0.0980;
              0.929 0.694 0.1250;
              0.494 0.184 0.556;
              0.466 0.674 0.188;
              0.301 0.7445 0.933;
              0.635 0.078 0.184];

    figure
    hold on
    title('First Story Column Base Acceptance Criteria: EW') 
    sct = scatter(frag.med_sa_ew,frag_io.prob_first,[],matlab_colors(1,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag_io.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First IO Column')
    sct = scatter(frag.med_sa_ew,frag_ls.prob_first,[],matlab_colors(2,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag_ls.frag_curve_ew_first,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','First LS Column')
    sct = scatter(frag.med_sa_ew,frag_cp.prob_first,[],matlab_colors(3,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag_cp.frag_curve_ew_first,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','First CP Column')
%     sct = scatter(frag.med_sa_ew,frag_euro.prob_first,[],color4(4,:),'filled');
%     set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     plot(x_points,frag_euro.frag_curve_ew_first,'color',color4(4,:),'DisplayName','First Euro Column')
    xlabel('Sa(T_1=1.14) (g)')
    ylabel('P[Exceedance]')
    xlim([0 1])
    fn_format_and_save_plot( plot_dir, 'Collapse Fragility - Acceptance Criteria - EW', 5 )
    
    % Save Collapse Fragilty Table
    
end

function [frag_curve_MLE, frag_curve_full_uncertainty] = fn_calc_frag_curve(x_points, med_sa, num_gms, num_collapse)

    [MLE_theta, MLE_beta] = fn_mle_pc(med_sa, num_gms, num_collapse);
    frag_curve_MLE = logncdf(x_points,log(MLE_theta),MLE_beta);
    frag_curve_full_uncertainty = logncdf(x_points,log(MLE_theta),0.6);
    
end

function [frag] = fn_exceedance_values(frag, i, num_occur, prcnt_occur)

    frag.num_first(i,1) = sum(num_occur >= 1);
    frag.num_10(i,1) = sum(prcnt_occur >= 0.1);
    frag.num_25(i,1) = sum(prcnt_occur >= 0.25);
    frag.num_50(i,1) = sum(prcnt_occur >= 0.5);
    frag.num_75(i,1) = sum(prcnt_occur >= 0.75);
    frag.num_100(i,1) = sum(prcnt_occur >= 1.0);
    frag.prob_first(i,1) = sum(num_occur >= 1) / length(num_occur);
    frag.prob_10(i,1) = sum(prcnt_occur >= 0.1) / length(prcnt_occur);
    frag.prob_25(i,1) = sum(prcnt_occur >= 0.25) / length(prcnt_occur);
    frag.prob_50(i,1) = sum(prcnt_occur >= 0.5) / length(prcnt_occur);
    frag.prob_75(i,1) = sum(prcnt_occur >= 0.75) / length(prcnt_occur);
    frag.prob_100(i,1) = sum(prcnt_occur >= 1.0) / length(prcnt_occur);
end

function [frag] = fn_multi_frag_curves(x_points, frag, med_sa_ew, num_gms)
    
    [frag.frag_curve_ew_first, frag.frag_curve_full_ew_first] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_first);
    [frag.frag_curve_ew_10, frag.frag_curve_full_ew_10] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_10);
    [frag.frag_curve_ew_25, frag.frag_curve_full_ew_25] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_25);
    [frag.frag_curve_ew_50, frag.frag_curve_full_ew_50] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_50);
    [frag.frag_curve_ew_75, frag.frag_curve_full_ew_75] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_75);
    [frag.frag_curve_ew_100, frag.frag_curve_full_ew_100] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_100);
    
end

function [ ] = fn_plot_frag_curves(x_points, frag, med_sa_ew, col_frag_curve, col_frag_curve_drift, icsb_motion, plot_name, plot_dir)    
% Import packages
import plotting_tools.fn_format_and_save_plot

%     color6 = [140,81,10;
%               216,179,101;
%               246,232,195;
%               199,234,229;
%               90,180,172;
%               1,102,94]/255;
          
%     color6 = [166 206 227;
%               31 120 180;
%               175 223 138;
%               51 160 44;
%               251 154 153;
%               227 26 28]/255;
%           
%           
                    
    matlab_colors = [0 0.4470 0.7410;
              0.85 0.325 0.0980;
              0.929 0.694 0.1250;
              0.494 0.184 0.556;
              0.466 0.674 0.188;
              0.301 0.7445 0.933;
              0.635 0.078 0.184];
              
    % Plot 1
    sct = scatter(med_sa_ew,frag.prob_first,[],matlab_colors(1,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Column Base')
    sct = scatter(med_sa_ew,frag.prob_10,[],matlab_colors(2,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Column Bases')
    sct = scatter(med_sa_ew,frag.prob_25,[],matlab_colors(3,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Column Bases')
    sct = scatter(med_sa_ew,frag.prob_50,[],matlab_colors(4,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Column Bases')
    sct = scatter(med_sa_ew,frag.prob_75,[],matlab_colors(5,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Column Bases')
    sct = scatter(med_sa_ew,frag.prob_100,[],matlab_colors(6,:),'filled');
    set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    plot(x_points,frag.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Column Bases')
    plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
    xlabel('Sa(T_1=1.14) (g)')
    
    ylabel('P[Exceedance]')
    xlim([0 1])
    fn_format_and_save_plot( plot_dir, plot_name, 5 )
    
    % Plot 2
    hold on
    plot(x_points,frag.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Column Base')
    plot(x_points,frag.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Column Bases')
    plot(x_points,frag.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Column Bases')
    plot(x_points,frag.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Column Bases')
    plot(x_points,frag.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Column Bases')
    plot(x_points,frag.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Column Bases')
    plot(x_points,col_frag_curve,'k','lineWidth',2,'DisplayName','Collapse')
%     plot(x_points,col_frag_curve_drift,':k','lineWidth',1.5,'DisplayName','Collapse from Drift')
    plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
    xlabel('Sa(T_1=1.14) (g)')
    
    ylabel('P[Exceedance]')
    xlim([0 1])
    fn_format_and_save_plot( plot_dir, [plot_name ' - plot2'], 5 )
    
end
