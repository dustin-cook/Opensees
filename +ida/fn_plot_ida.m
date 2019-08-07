function [ ] = fn_plot_ida(analysis, model, IDA_scale_factors, gm_set_table, max_dir_spectra, ida_results, SSF_ew, SSF_ns)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.fn_format_and_save_plot

% Defined fixed parames
params = {'b','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};

% Collect IDA data
plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'IDA Plots'];
if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end
id = 0;
id_missing = 0;
for i = 1:length(IDA_scale_factors)
    for gms = 1:height(gm_set_table)   
        % Load data
        outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(IDA_scale_factors(i)) '/' 'GM_' num2str(gm_set_table.set_id(gms)) '_' num2str(gm_set_table.pair(gms))];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
        if exist(outputs_file,'file') && exist(hinge_file,'file')
            id = id + 1;
            load(outputs_file)
            load(hinge_file)
            ida.id(id,1) = id;
            ida.scale(id,1) = IDA_scale_factors(i);
            ida.gm_name{id,1} = gm_set_table.eq_name{gms};
            
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
            
            % Collapse metrics
            ida.collapse(id,1) = summary.collapse;
            if isfield(summary,'collapse_direction')
                ida.collapse_direction{id,1} = summary.collapse_direction;
            end
            if isfield(summary,'collaspe_mech')
                ida.collapse_mech{id,1} = summary.collaspe_mech;
            end
            
            % Get element group filters
            col_base_filter = hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column');
            first_story_col_filter = hinge.story == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column');
            first_story_beam_filter = hinge.story == 1 & strcmp(hinge.direction,'primary')  & strcmp(hinge.ele_type,'beam');
            col_hinges_base = hinge(col_base_filter,:);
            col_hinges_story_1 = hinge(first_story_col_filter,:);
            bm_hinges_story_1 = hinge(first_story_beam_filter,:);

            % For Each accetance criteria listed above
            for p = 1:length(params)
                [ ida.(['num_' params{p}])(id,1), ida.(['percent_' params{p}])(id,1), ida.(['num_' params{p} '_15'])(id,1) ] = fn_collect_ida_data([params{p} '_ratio'], col_hinges_base, summary.collapse);
                [ ida.(['num_' params{p} '_cols_1'])(id,1), ida.(['percent_' params{p} '_cols_1'])(id,1), ida.(['num_' params{p} '_15' '_cols_1'])(id,1) ] = fn_collect_ida_data([params{p} '_ratio'], col_hinges_story_1, summary.collapse);
                [ ida.(['num_' params{p} '_bms_1'])(id,1), ida.(['percent_' params{p} '_bms_1'])(id,1), ida.(['num_' params{p} '_15' '_bms_1'])(id,1) ] = fn_collect_ida_data([params{p} '_ratio'], bm_hinges_story_1, summary.collapse);
            end

            % Shear ASCE 41
            col_hinges_base.shear_ratio = col_hinges_base.shear_demand ./ col_hinges_base.asce41_shear_capacity;
            [ ida.num_shear_asce41(id,1), ida.percent_shear_asce41(id,1), ~ ] = fn_collect_ida_data('shear_ratio', col_hinges_base, summary.collapse);

            % Shear Eurocode
            ida.shear_asce41_v_euro(id,1) = sum((col_hinges_base.euro_V_NC ./ col_hinges_base.asce41_shear_capacity) > 1);

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
if exist('missing_ida','var')
    writetable(struct2table(missing_ida),[plot_dir filesep 'idas_missing.csv'])
end
writetable(failed_convergence,[plot_dir filesep 'idas_failed_convergence.csv'])

%% Calculate Median Drifts and Accels
for i = 1:length(IDA_scale_factors)
    %% X Direction
    % Lognormal Mean Spectral Accelration
    frag.p695_sa_ew(i,1) = exp(mean(log(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)))));
    frag.p_15_sa_ew(i,1) = prctile(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)),15);
    frag.p_85_sa_ew(i,1) = prctile(ida_table.sa_x(ida_table.scale == IDA_scale_factors(i)),85);

    % Max Direction Spectral Acceleration
    [~,idx_ew] = min(abs(max_dir_spectra.period - ida_results.period(1)));
    frag.asce41_sa_ew(i,1) = max_dir_spectra.sa(idx_ew)*IDA_scale_factors(i);

    % Drift
    frag.mean_idr_ew(i,1) = mean(ida_table.drift_x(ida_table.scale == IDA_scale_factors(i)));
    
    %% Z Direction
    if analysis.run_z_motion
        % Lognormal Mean Spectral Accelration
        frag.p695_sa_ns(i,1) = exp(mean(log(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)))));
        frag.p_15_sa_ns(i,1) = prctile(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)),15);
        frag.p_85_sa_ns(i,1) = prctile(ida_table.sa_z(ida_table.scale == IDA_scale_factors(i)),85);

        % Max Direction Spectral Acceleration
        [~,idx_ns] = min(abs(max_dir_spectra.period - ida_results.period(2)));
        frag.asce41_sa_ns(i,1) = max_dir_spectra.sa(idx_ns)*IDA_scale_factors(i);

        % Drift
        frag.mean_idr_ns(i,1) = exp(mean(log(ida_table.drift_z(ida_table.scale == IDA_scale_factors(i)))));
    end
end

%% Plot IDA curves
fprintf('Saving IDA Summary Data and Figures to Directory: %s',plot_dir)
% Plot X direction IDA curves
hold on
title('Plan EW (x)')
for gms = 1:height(gm_set_table)
    ida_plt = plot(ida_table.drift_x(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),ida_table.sa_x(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),'color',[0.75 0.75 0.75]);
    set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
end
plot(frag.mean_idr_ew,frag.p695_sa_ew,'b','lineWidth',1.5,'DisplayName','Mean Drift')
plot(frag.mean_idr_ew,frag.p_15_sa_ew,'--b','lineWidth',1.5,'DisplayName','15th Percentile')
plot(frag.mean_idr_ew,frag.p_85_sa_ew,'--b','lineWidth',1.5,'DisplayName','85th Percentile')
plot([0,0.08],[ida_results.spectra(1),ida_results.spectra(1)],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
xlim([0 0.08])
xlabel('Max Drift')
ylabel('Sa(T_1=1.14s) (g)')
fn_format_and_save_plot( plot_dir, 'IDA Plot EW Frame Direction', 3 )

% Plot Z direction IDA curves
if analysis.run_z_motion
    hold on
    title('Plan NS (y)')
    for gms = 1:height(gm_set_table)
        ida_plt = plot(ida_table.drift_z(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),ida_table.sa_z(strcmp(ida_table.gm_name,gm_set_table.eq_name{gms})),'color',[0.75 0.75 0.75]);
        set(get(get(ida_plt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    plot(frag.mean_idr_ns,frag.p695_sa_ns,'b','lineWidth',1.5,'DisplayName','Mean Drift')
    plot(frag.mean_idr_ns,frag.p_15_sa_ns,'--b','lineWidth',1.5,'DisplayName','15th Percentile')
    plot(frag.mean_idr_ns,frag.p_85_sa_ns,'--b','lineWidth',1.5,'DisplayName','85th Percentile')
    plot([0,0.08],[ida_results.spectra(2),ida_results.spectra(2)],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
    xlim([0 0.08])
    xlabel('Max Drift')
    ylabel('Sa(T_1=0.44s) (g)')
    fn_format_and_save_plot( plot_dir, 'IDA Plot NS Wall Direction', 3 )
end

%% Collect frag curve info
% Collect Stripe info
frag_b = [];
frag_io = [];
frag_ls = [];
frag_cp = [];
frag_euro_NC = [];
frag_euro_SD = [];
frag_euro_DL = [];
frag_b_drift = [];
frag_shear_asce41 = [];
for i = 1:length(IDA_scale_factors)
    stripe = ida_table(ida_table.scale == IDA_scale_factors(i),:);

    % Regular Collapse
    frag.num_gms(i,1) = height(stripe);
    frag.num_gms_drift(i,1) = sum(stripe.collapse ~= 3);
    frag.num_gms_convergence(i,1) = sum(stripe.collapse ~= 1);
    frag.num_collapse_drift(i,1) = sum(stripe.collapse == 1);
    frag.num_collapse_convergence(i,1) = sum(stripe.collapse == 3);
    frag.num_collapse_tot(i,1) = frag.num_collapse_drift(i,1) + frag.num_collapse_convergence(i,1);
    frag.num_b_15(i,1) = sum(stripe.num_b_15 > 0);
    frag.num_unacceptable_response(i,1) = sum(stripe.collapse == 1 | stripe.collapse == 3 | stripe.num_b_15 > 0);
    frag.prob_collapse_drift(i,1) = frag.num_collapse_drift(i,1) / frag.num_gms(i,1);
    frag.prob_collapse_convergence(i,1) = frag.num_collapse_convergence(i,1) / frag.num_gms(i,1);
    frag.prob_collapse_tot(i,1) = frag.num_collapse_tot(i,1) / frag.num_gms(i,1);
    frag.prob_b_15(i,1) = frag.num_b_15(i,1) / frag.num_gms(i,1);
    frag.prob_unacceptable_response(i,1) = frag.num_unacceptable_response(i,1) / frag.num_gms(i,1);

    % For Each accetance criteria listed above
    for p = 1:length(params)
        [new_row] = fn_exceedance_values(stripe.(['num_' params{p}]), stripe.(['percent_' params{p}]));
        frag.(params{p})(i,:) = new_row;
    end

    % All fist story columns
    [new_row] = fn_exceedance_values(stripe.(['num_b_cols_1']), stripe.(['percent_b_cols_1']));
    frag.b_cols_1(i,:) = new_row;
            
    % Shear Capacity
    [new_row] = fn_exceedance_values(stripe.num_shear_asce41, stripe.percent_shear_asce41);
    frag.shear_asce41(i,:) = new_row;
end

%% Create Fragility Curves based on Baker MLE
x_points = 0:0.01:3;
% X direction Curves
[col_frag_curve_ew_drift, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ew, frag.num_gms_drift, frag.num_collapse_drift);
[col_frag_curve_ew_convergence, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ew, frag.num_gms_convergence, frag.num_collapse_convergence);
[col_frag_curve_ew_b_15, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ew, frag.num_gms, frag.num_b_15);
[col_frag_curve_ew_tot, col_frag_curve_full_ew_tot] = fn_calc_frag_curve(x_points, frag.asce41_sa_ew, frag.num_gms, frag.num_collapse_tot);
[col_frag_curve_ew_unacceptable_response, col_frag_curve_full_ew_unacceptable_response] = fn_calc_frag_curve(x_points, frag.asce41_sa_ew, frag.num_gms, frag.num_unacceptable_response);
for p = 1:length(params)
    [frag_curves.(params{p})] = fn_multi_frag_curves(x_points, frag.(params{p}), frag.asce41_sa_ew, frag.num_gms);
end
[frag_curves.b_cols_1] = fn_multi_frag_curves(x_points, frag.b_cols_1, frag.asce41_sa_ew, frag.num_gms);
[frag_curves.shear_asce41] = fn_multi_frag_curves(x_points, frag.shear_asce41, frag.asce41_sa_ew, frag.num_gms);
[frag_curves.(params{p})] = fn_multi_frag_curves(x_points, frag.(params{p}), frag.asce41_sa_ew, frag.num_gms);

% Z direction curves
if analysis.run_z_motion
    [col_frag_curve_ns_drift, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ns, frag.num_gms_drift, frag.num_collapse_drift);
    [col_frag_curve_ns_convergence, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ns, frag.num_gms_convergence, frag.num_collapse_convergence);
    [col_frag_curve_ns_b_15, ~] = fn_calc_frag_curve(x_points, frag.asce41_sa_ns, frag.num_gms, frag.num_b_15);
    [col_frag_curve_ns_tot, col_frag_curve_full_ns_tot] = fn_calc_frag_curve(x_points, frag.asce41_sa_ns, frag.num_gms, frag.num_collapse_tot);
    [col_frag_curve_ns_unacceptable_response, col_frag_curve_full_ns_unacceptable_response] = fn_calc_frag_curve(x_points, frag.asce41_sa_ns, frag.num_gms, frag.num_unacceptable_response);
end
%% Calculate Post Fragulity Curve P695 factors
if analysis.run_z_motion
    factor_3D = 1.2;
else
    factor_3D = 1.0;
end
% X direction
% Collapse Median Sa
[~, idx_med_ew] = min(abs(col_frag_curve_ew_tot - 0.5));
ida_results.sa_med_col(1,1) = x_points(idx_med_ew);

% Collapse Margin Ratio
ida_results.cmr(1,1) = ida_results.sa_med_col(1) / ida_results.mce(1);

% Adjust for SSF and 3D
ida_results.acmr(1,1) = factor_3D*SSF_ew*ida_results.cmr(1);
median_adjustment(1) = ida_results.sa_med_col(1)*(factor_3D*SSF_ew - 1);


% Create Fragility curves based on P695
x_points_adjused_x = x_points + median_adjustment(1);

% Find P[C|MCE]
[~, idx_ew] = min(abs(x_points - ida_results.mce(1)));
ida_results.p_col_mce(1,1) = col_frag_curve_ew_tot(idx_ew);
[~, idx_ew] = min(abs(x_points_adjused_x - ida_results.mce(1)));
ida_results.p_col_mce_adjust(1,1) = col_frag_curve_full_ew_tot(idx_ew);
[~, idx_UR] = min(abs(x_points - ida_results.mce(1)));
p_UR_mce(1,1) = col_frag_curve_ew_unacceptable_response(idx_UR);
[~, idx_UR] = min(abs(x_points_adjused_x - ida_results.mce(1)));
p_UR_mce_adjusted(1,1) = col_frag_curve_full_ew_unacceptable_response(idx_UR);

% Z direction
if analysis.run_z_motion
    % Collapse Median Sa
    [~, idx_med_ns] = min(abs(col_frag_curve_ns_tot - 0.5));
    ida_results.sa_med_col(2,1) = x_points(idx_med_ns);

    % Collapse Margin Ratio
    ida_results.cmr(2,1) = ida_results.sa_med_col(2) / ida_results.mce(2);

    % Adjust for SSF and 3D
    ida_results.acmr(2,1) = factor_3D*SSF_ns*ida_results.cmr(2);
    median_adjustment(2) = ida_results.sa_med_col(2)*(factor_3D*SSF_ns - 1);

    % Create Fragility curves based on P695
    x_points_adjused_z = x_points + median_adjustment(2);

    % Find P[C|MCE]
    [~, idx_ns] = min(abs(x_points - ida_results.mce(2)));
    ida_results.p_col_mce(2,1) = col_frag_curve_ns_tot(idx_ns);
    [~, idx_ns] = min(abs(x_points_adjused_z - ida_results.mce(2)));
    ida_results.p_col_mce_adjust(2,1) = col_frag_curve_full_ns_tot(idx_ns);
    [~, idx_UR] = min(abs(x_points - ida_results.mce(2)));
    p_UR_mce(2,1) = col_frag_curve_ns_unacceptable_response(idx_UR);
    [~, idx_UR] = min(abs(x_points_adjused_x - ida_results.mce(2)));
    p_UR_mce_adjusted(2,1) = col_frag_curve_full_ns_unacceptable_response(idx_UR);
end

% Save IDA results as table
ida_results_table = struct2table(ida_results);
writetable(ida_results_table,[plot_dir filesep 'ida_results.csv'])

%% Plot Frag Curves    
% Total Collapse with discrete scatter
figure
hold on
title('P-695 Collapse Fragility') 
sct = scatter(frag.asce41_sa_ew,frag.prob_collapse_tot,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_tot,'b','DisplayName','EW Frame Fragility')
plot([ida_results.spectra(1),ida_results.spectra(1)],[0,1],'--b','lineWidth',1.0,'DisplayName','ICSB EW Spectra')
% sct = scatter(frag.asce41_sa_ns,frag.prob_collapse_tot,'k','filled');
% set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
% plot(x_points,col_frag_curve_ns_tot,'k','DisplayName','NS Wall Fragility')
% plot([ida_results.spectra(2),ida_results.spectra(2)],[0,1],'--k','lineWidth',1.0,'DisplayName','ICSB NS Spectra')
xlabel('Sa(T_1) (g)')
ylabel('P[Collapse]')
fn_format_and_save_plot( plot_dir, 'Collapse Fragility', 6 )

% Total Collapse with 695 adjustments
figure
hold on
title('P-695 Collapse Fragility') 
plot(x_points,col_frag_curve_ew_tot,'b','DisplayName','EW Frame Fragility')
plot(x_points_adjused_x,col_frag_curve_ew_tot,'--b','DisplayName','EW SSF Adjustment')
plot(x_points_adjused_x,col_frag_curve_full_ew_tot,':b','DisplayName','EW Full Uncertainty')
% plot(x_points,col_frag_curve_ns_tot,'k','DisplayName','NS Wall Fragility')
% plot(x_points_adjused_z,col_frag_curve_ns_tot,'--k','DisplayName','NS SSF Adjustment')
% plot(x_points_adjused_z,col_frag_curve_full_ns_tot,':k','DisplayName','NS Full Uncertainty')
xlabel('Sa(T_1) (g)')
ylabel('P[Collapse]')
fn_format_and_save_plot( plot_dir, 'Collapse Fragility with 695 Uncertainty', 6 )

% Discrete Collapse Breakdown - EW
figure
hold on
title('P-695 Collapse Fragility: EW Frame') 
plot(frag.asce41_sa_ew,frag.prob_collapse_tot,'k','DisplayName','Total');
plot(frag.asce41_sa_ew,frag.prob_collapse_drift,'b','DisplayName','6% Drift');
plot(frag.asce41_sa_ew,frag.prob_collapse_convergence,'r','DisplayName','Nonconvergence');
xlabel('Sa(T_1=1.14) (g)')
ylabel('P[Collapse]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, 'Discrete Collapse Fragility Breakdown - EW', 6 )

% Continous Collapse Breakdown - EW
figure
hold on
title('P-695 Collapse Fragility: EW Frame') 
sct = scatter(frag.asce41_sa_ew,frag.prob_collapse_tot,'k','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_tot,'k','DisplayName','Total')
sct = scatter(frag.asce41_sa_ew,frag.num_collapse_drift./frag.num_gms_drift,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_drift,'b','DisplayName','6% Drift')
sct = scatter(frag.asce41_sa_ew,frag.num_collapse_convergence./frag.num_gms_convergence,'r','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_convergence,'r','DisplayName','Nonconvergence')
xlabel('Sa(T_1) (g)')
ylabel('P[Collapse]')
xlim([0,1.5])
fn_format_and_save_plot( plot_dir, 'Continuous Collapse Fragility Breakdown - EW', 6 )

% ASCE 41 Unacceptable Response v Total Collapse with 695 adjustments
figure
hold on
title('ASCE 41 Unacceptable Response') 
plot(x_points,col_frag_curve_ew_unacceptable_response,'b','DisplayName','EW Unacceptable Response')
plot(x_points,col_frag_curve_ew_tot,'--b','DisplayName','EW Collapse')
plot(x_points_adjused_x,col_frag_curve_full_ew_tot,':b','DisplayName','EW Adjusted Collapse')
% plot(x_points,col_frag_curve_ns_unacceptable_response,'k','DisplayName','NS Unacceptable Response')
% plot(x_points,col_frag_curve_ns_tot,'--k','DisplayName','NS Collapse')
% plot(x_points_adjused_z,col_frag_curve_full_ns_tot,':k','DisplayName','NS Adjusted Collapse')
xlabel('Sa(T_1) (g)')
ylabel('P[Exceedance]')
fn_format_and_save_plot( plot_dir, 'Unacceptable Response Fragility', 6 )

% Discrete Unacceptable response breakdown - EW
figure
hold on
title('ASCE 41 Unacceptable Response: EW Frame') 
plot(frag.asce41_sa_ew,frag.prob_unacceptable_response,'k','DisplayName','Total');
plot(frag.asce41_sa_ew,frag.prob_b_15,'m','DisplayName','Element Rotation > 1.5*b');
plot(frag.asce41_sa_ew,frag.prob_collapse_drift,'b','DisplayName','6% Drift');
plot(frag.asce41_sa_ew,frag.prob_collapse_convergence,'r','DisplayName','Nonconvergence');
xlabel('Sa(T_1=1.14) (g)')
ylabel('P[Unacceptable Response]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, 'Discrete Unacceptable Response Breakdown - EW', 6 )

% Continous Unacceptable Response Breakdown - EW
figure
hold on
title('ASCE 41 Unacceptable Response: EW Frame') 
sct = scatter(frag.asce41_sa_ew,frag.prob_unacceptable_response,'k','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_unacceptable_response,'k','DisplayName','Total')
sct = scatter(frag.asce41_sa_ew,frag.prob_b_15,'m','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_b_15,'m','DisplayName','Element Rotation > 1.5*b')
sct = scatter(frag.asce41_sa_ew,frag.num_collapse_drift./frag.num_gms_drift,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_drift,'b','DisplayName','6% Drift')
sct = scatter(frag.asce41_sa_ew,frag.num_collapse_convergence./frag.num_gms_convergence,'r','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,col_frag_curve_ew_convergence,'r','DisplayName','Nonconvergence')
xlabel('Sa(T_1) (g)')
ylabel('P[Unacceptable Response]')
xlim([0,1.5])
fn_format_and_save_plot( plot_dir, 'Continuous Unacceptable Response Breakdown - EW', 6 )

% Plot Each acceptance Criteria
for p = 1:length(params)
    plot_title = ['First Story Column Base ' params{p} ' Fragility: EW'];
    plot_name = ['Collapse Fragility - ' params{p} ' - EW'];
    fn_plot_frag_curves(x_points, frag.(params{p}), frag_curves.(params{p}), frag.asce41_sa_ew, col_frag_curve_ew_tot, ida_results.spectra(1), plot_name, plot_dir, plot_title)
end

% Columns both sides - first story mechanism
plot_title = 'First Story Columns b Fragility: EW';
plot_name = 'Collapse Fragility - first story b - EW';
fn_plot_frag_curves(x_points, frag.b_cols_1, frag_curves.b_cols_1, frag.asce41_sa_ew, col_frag_curve_ew_tot, ida_results.spectra(1), plot_name, plot_dir, plot_title)

% Column Shear ASCE 41
plot_title = 'First Story Column Base ASCE 41 Shear Capacity Fragility: EW';
plot_name = 'Collapse Fragility - asce41 shear - EW';
fn_plot_frag_curves(x_points, frag.shear_asce41, frag_curves.shear_asce41, frag.asce41_sa_ew, col_frag_curve_ew_tot, ida_results.spectra(1), plot_name, plot_dir, plot_title)

% Combined acceptance plots  
matlab_colors = [0 0.4470 0.7410;
          0.85 0.325 0.0980;
          0.929 0.694 0.1250;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.7445 0.933;
          0.635 0.078 0.184];

% First Acceptance criteria met
figure
hold on
title('First Story Column Base Acceptance Criteria: EW') 
% ASCE 41 frags
plot(x_points,frag_curves.io.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First IO')
plot(x_points,frag_curves.ls.frag_curve_ew_first,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','First LS')
plot(x_points,frag_curves.cp.frag_curve_ew_first,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','First CP')
% Eurcode Frags
plot(x_points,frag_curves.euro_th_DL.frag_curve_ew_first,'--','color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Euro DL')
plot(x_points,frag_curves.euro_th_SD.frag_curve_ew_first,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','First Euro SD')
plot(x_points,frag_curves.euro_th_NC.frag_curve_ew_first,'--','color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','First Euro NC')
xlabel('Sa(T_1=1.14) (g)')
ylabel('P[Exceedance]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility - First Acceptance Criteria - EW', 6 )

% 50% Acceptance criteria met
figure
hold on
title('50% Story Column Base Acceptance Criteria: EW') 
% ASCE 41 frags
plot(x_points,frag_curves.io.frag_curve_ew_50,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','50% IO')
plot(x_points,frag_curves.ls.frag_curve_ew_50,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','50% LS')
plot(x_points,frag_curves.cp.frag_curve_ew_50,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','50% CP')
% Eurcode Frags
plot(x_points,frag_curves.euro_th_DL.frag_curve_ew_50,'--','color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','50% Euro DL')
plot(x_points,frag_curves.euro_th_SD.frag_curve_ew_50,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','50% Euro SD')
plot(x_points,frag_curves.euro_th_NC.frag_curve_ew_50,'--','color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','50% Euro NC')
xlabel('Sa(T_1=1.14) (g)')
ylabel('P[Exceedance]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility - 50% Acceptance Criteria - EW', 6 )
end

function [frag_curve_MLE, frag_curve_full_uncertainty] = fn_calc_frag_curve(x_points, med_sa, num_gms, num_collapse)
import ida.fn_mle_pc
[MLE_theta, MLE_beta] = fn_mle_pc(med_sa, num_gms, num_collapse);
frag_curve_MLE = logncdf(x_points,log(MLE_theta),MLE_beta);
frag_curve_full_uncertainty = logncdf(x_points,log(MLE_theta),0.6);
end

function [frag_table] = fn_exceedance_values(num_occur, prcnt_occur)

frag.num_first = sum(num_occur >= 1);
frag.num_10 = sum(prcnt_occur >= 0.1);
frag.num_25 = sum(prcnt_occur >= 0.25);
frag.num_50 = sum(prcnt_occur >= 0.5);
frag.num_75 = sum(prcnt_occur >= 0.75);
frag.num_100 = sum(prcnt_occur >= 1.0);
frag.prob_first = sum(num_occur >= 1) / length(num_occur);
frag.prob_10 = sum(prcnt_occur >= 0.1) / length(prcnt_occur);
frag.prob_25 = sum(prcnt_occur >= 0.25) / length(prcnt_occur);
frag.prob_50 = sum(prcnt_occur >= 0.5) / length(prcnt_occur);
frag.prob_75 = sum(prcnt_occur >= 0.75) / length(prcnt_occur);
frag.prob_100 = sum(prcnt_occur >= 1.0) / length(prcnt_occur);
frag_table = struct2table(frag);
end

function [frag_curves] = fn_multi_frag_curves(x_points, frag, med_sa_ew, num_gms)
    
[frag_curves.frag_curve_ew_first, frag_curves.frag_curve_full_ew_first] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_first);
[frag_curves.frag_curve_ew_10, frag_curves.frag_curve_full_ew_10] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_10);
[frag_curves.frag_curve_ew_25, frag_curves.frag_curve_full_ew_25] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_25);
[frag_curves.frag_curve_ew_50, frag_curves.frag_curve_full_ew_50] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_50);
[frag_curves.frag_curve_ew_75, frag_curves.frag_curve_full_ew_75] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_75);
[frag_curves.frag_curve_ew_100, frag_curves.frag_curve_full_ew_100] = fn_calc_frag_curve(x_points, med_sa_ew, num_gms, frag.num_100);
end

function [ ] = fn_plot_frag_curves(x_points, frag, frag_curves, med_sa_ew, col_frag_curve, icsb_motion, plot_name, plot_dir, plot_title)    
% Import packages
import plotting_tools.fn_format_and_save_plot
                
matlab_colors = [0 0.4470 0.7410;
          0.85 0.325 0.0980;
          0.929 0.694 0.1250;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.7445 0.933;
          0.635 0.078 0.184];

% Plot 1
hold on
title(plot_title)
sct = scatter(med_sa_ew,frag.prob_first,[],matlab_colors(1,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Column Base')
sct = scatter(med_sa_ew,frag.prob_10,[],matlab_colors(2,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Column Bases')
sct = scatter(med_sa_ew,frag.prob_25,[],matlab_colors(3,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Column Bases')
sct = scatter(med_sa_ew,frag.prob_50,[],matlab_colors(4,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Column Bases')
sct = scatter(med_sa_ew,frag.prob_75,[],matlab_colors(5,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Column Bases')
sct = scatter(med_sa_ew,frag.prob_100,[],matlab_colors(6,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Column Bases')
plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
xlabel('Sa(T_1=1.14) (g)')

ylabel('P[Exceedance]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, plot_name, 6 )

% Plot 2
hold on
title(plot_title)
plot(x_points,frag_curves.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Element')
plot(x_points,frag_curves.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Elements')
plot(x_points,frag_curves.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Elements')
plot(x_points,frag_curves.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Elements')
plot(x_points,frag_curves.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Elements')
plot(x_points,frag_curves.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Elements')
plot(x_points,col_frag_curve,'k','lineWidth',2,'DisplayName','Collapse')
plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
xlabel('Sa(T_1=1.14) (g)')

ylabel('P[Exceedance]')
xlim([0 1.5])
fn_format_and_save_plot( plot_dir, [plot_name ' - plot2'], 6 )
end

function [ num_cols, percent_cols, num_cols_15 ] = fn_collect_ida_data(var_name, col_hinges, collapse)

col_ratios = col_hinges.(var_name);
if collapse == 3 || collapse == 1
    num_cols = length(col_ratios);
else
    num_cols = sum(col_ratios >= 1.0);
end
percent_cols = num_cols / length(col_ratios);
num_cols_15 = sum(col_ratios >= 1.5);
end