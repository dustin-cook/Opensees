function [ ] = fn_plot_ida(analysis, model, gm_set_table, ida_results, SSF_ew)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};

% Define Directiries
read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end

% Load data
ida_table = readtable([read_dir filesep 'ida_table.csv']);
load([read_dir filesep 'frag_data.mat'])
load([read_dir filesep 'frag_curves.mat'])

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

%% Calculate Post Fragulity Curve P695 factors
if analysis.run_z_motion
    factor_3D = 1.11;
else
    factor_3D = 1.0;
end
SSF = SSF_ew; % Assume EW SSF

% Set Ida summary table 
ida_summary_table.direction = 'EW';
ida_summary_table.period = ida_results.period(1);
ida_summary_table.spectra = ida_results.spectra(1);
ida_summary_table.mce = ida_results.mce(1);

% Collapse Median Sa
[~, idx_med] = min(abs(frag_curves.collapse_tot - 0.5));
ida_summary_table.sa_med_col = x_points(idx_med);
[~, idx_med] = min(abs(frag_curves.unacceptable_response - 0.5));
ida_summary_table.sa_med_UR = x_points(idx_med);
[~, idx_med] = min(abs(frag_curves.cp_cols_walls_1.frag_curve_ew_first - 0.5));
ida_summary_table.sa_med_CP = x_points(idx_med);

% Collapse Margin Ratio
ida_summary_table.cmr = ida_summary_table.sa_med_col / ida_summary_table.mce;
ida_summary_table.cmr_UR = ida_summary_table.sa_med_UR / ida_summary_table.mce;
ida_summary_table.cmr_CP = ida_summary_table.sa_med_CP / ida_summary_table.mce;

% Adjust for SSF and 3D
ida_summary_table.acmr = factor_3D*SSF*ida_summary_table.cmr;
ida_summary_table.acmr_UR = factor_3D*SSF*ida_summary_table.cmr_UR;
ida_summary_table.acmr_CP = factor_3D*SSF*ida_summary_table.cmr_CP;
median_adjustment = ida_summary_table.sa_med_col*(factor_3D*SSF - 1);

% Create Fragility curves based on P695
x_points_adjused_x = x_points + median_adjustment(1);

% Find P[C|MCE]
[~, idx] = min(abs(x_points - ida_summary_table.mce));
[~, idx_adjust] = min(abs(x_points_adjused_x - ida_summary_table.mce));
[~, idx_ICSB] = min(abs(x_points - ida_summary_table.spectra));

ida_summary_table.p_col_mce = frag_curves.collapse_tot(idx);
ida_summary_table.p_col_mce_adjust = frag_curves.collapse_full_tot(idx_adjust);
ida_summary_table.p_UR_mce = frag_curves.unacceptable_response(idx);
ida_summary_table.p_UR_mce_adjusted = frag_curves.full_unacceptable_response(idx_adjust);
ida_summary_table.p_CP_mce = frag_curves.cp_cols_walls_1.frag_curve_ew_first(idx);
ida_summary_table.p_C_ICSM = frag_curves.collapse_tot(idx_ICSB);

% Save IDA results as table
ida_results_table = struct2table(ida_summary_table);
writetable(ida_results_table,[plot_dir filesep 'ida_results.csv'])

%% Plot Frag Curves
% color scheme 
matlab_colors = [0 0.4470 0.7410;
          0.85 0.325 0.0980;
          0.929 0.694 0.1250;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.7445 0.933;
          0.635 0.078 0.184];
      
% Total Collapse with discrete scatter
figure
hold on
% title('P-695 Collapse Fragility') 
sct = scatter(sa_points,frag.num_collapse_tot./frag.num_gms,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_tot,'b','DisplayName','Collapse Fragility')
% plot([ida_results.spectra(1),ida_results.spectra(1)],[0,1],'--b','lineWidth',1.0,'DisplayName','ICSB EW Spectra')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Collapse]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility', 6 )

% Total Collapse with 695 adjustments
figure
hold on
% title('P-695 Collapse Fragility') 
plot(x_points,frag_curves.collapse_tot,'b','lineWidth',1,'DisplayName','Raw Collapse')
plot(x_points_adjused_x,frag_curves.collapse_tot,'--b','lineWidth',1,'DisplayName','Adjustment')
plot(x_points_adjused_x,frag_curves.collapse_full_tot,':b','lineWidth',1.25,'DisplayName','With Uncertainty')
plot([ida_results.mce(1),ida_results.mce(1)],[0,1],'--k','lineWidth',0.5,'DisplayName','MCE_R')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Collapse]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility with 695 Uncertainty', 6 )

% Total Collapse with breakdown in each direction
figure
hold on
% title('P-695 Collapse Fragility') 
sct = scatter(sa_points,frag.num_collapse_tot ./ frag.num_gms,'k','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_tot,'k','DisplayName','Total Collapse Fragility')
sct = scatter(sa_points,frag.num_collapse_ew ./ frag.num_gms_ew,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_ew,'b','DisplayName','EW Frame Fragility')
sct = scatter(sa_points,frag.num_collapse_ns ./ frag.num_gms_ns,'r','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_ns,'r','DisplayName','NS Wall Fragility')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Collapse]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility Both Directions', 6 )

% Continous Collapse Breakdown 
figure
hold on
% title('P-695 Collapse Fragility') 
sct = scatter(sa_points,frag.num_collapse_tot ./ frag.num_gms,'k','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_tot,'k','DisplayName','Total')
sct = scatter(sa_points,frag.num_collapse_drift./frag.num_gms_drift,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_drift,'b','DisplayName','6% Drift')
sct = scatter(sa_points,frag.num_collapse_convergence./frag.num_gms_convergence,'r','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_convergence,'r','DisplayName','Nonconvergence')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Collapse]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Collapse Fragility Breakdown', 6 )

% ASCE 41 Unacceptable Response v Total Collapse with 695 adjustments
figure
hold on
% title('ASCE 41 Unacceptable Response') 
plot(x_points,frag_curves.unacceptable_response,'b','DisplayName','Unacceptable Response')
plot(x_points,frag_curves.collapse_tot,'--b','DisplayName','Collapse Fragility')
plot([ida_results.mce(1),ida_results.mce(1)],[0,1],'--k','lineWidth',0.5,'DisplayName','MCE_R')
% plot(x_points_adjused_x,col_frag_curve_full_tot,':b','DisplayName','Adjusted Collapse')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Unacceptable Response Fragility', 6 )

% Unacceptable Response Breakdown
figure
hold on
% title('ASCE 41 Unacceptable Response') 
sct = scatter(sa_points,frag.num_unacceptable_response ./ frag.num_gms,'k','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.unacceptable_response,'k','DisplayName','Total')
sct = scatter(sa_points,frag.num_accept_15 ./ frag.num_gms,'m','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_accept_15,'m','DisplayName','> 1.5 Element Failure')
sct = scatter(sa_points,frag.num_collapse_drift./frag.num_gms_drift,'b','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_drift,'b','DisplayName','6% Drift')
sct = scatter(sa_points,frag.num_collapse_convergence./frag.num_gms_convergence,'r','filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.collapse_convergence,'r','DisplayName','Nonconvergence')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Unacceptable Response]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'Unacceptable Response Breakdown', 6 )

% Plot Each acceptance Criteria
for p = 1:length(params)
    plot_title = ['First Story Column ' params{p} ' Fragilities'];
    plot_name = ['First Story Column ' params{p} ' Fragilities'];
    fn_plot_frag_curves(x_points, frag_curves.([params{p} '_cols_1']), frag_curves.collapse_ew, ida_results.spectra(1), plot_name, plot_dir, plot_title, [0,ceil(max(sa_points))], matlab_colors, frag.([params{p} '_cols_1']), frag.num_gms_ew, sa_points )

    plot_title = ['First Story Wall ' params{p} ' Fragilities'];
    plot_name = ['First Story Wall ' params{p} ' Fragilities'];
    fn_plot_frag_curves(x_points, frag_curves.([params{p} '_walls_1']), frag_curves.collapse_ns, ida_results.spectra(1), plot_name, plot_dir, plot_title, [0,ceil(max(sa_points))], matlab_colors, frag.([params{p} '_walls_1']), frag.num_gms_ns, sa_points )
    
    plot_title = ['First Story Mechanism ' params{p} ' Fragilities'];
    plot_name = ['First Story Mechanism ' params{p} ' Fragilities'];
    fn_plot_frag_curves(x_points, frag_curves.([params{p} '_cols_walls_1']), frag_curves.collapse_tot, ida_results.spectra(1), plot_name, plot_dir, plot_title, [0,ceil(max(sa_points))], matlab_colors, frag.([params{p} '_cols_walls_1']), frag.num_gms, sa_points )
end

% First Acceptance criteria met
figure
hold on
% title('First Acceptance Criteria Triggered') 
% ASCE 41 frags
% plot(x_points,frag_curves.io_cols_1.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 IO')
% plot(x_points,frag_curves.ls_cols_1.frag_curve_ew_first,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','ASCE 41 LS')
plot(x_points,frag_curves.cp_cols_1.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 CP')
% Eurcode Frags
% plot(x_points,frag_curves.euro_th_DL_cols_1.frag_curve_ew_first,'--','color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','EuroCode DL')
% plot(x_points,frag_curves.euro_th_SD_cols_1.frag_curve_ew_first,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode SD')
plot(x_points,frag_curves.euro_th_NC_cols_1.frag_curve_ew_first,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode NC')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'EuroCode - First Acceptance Criteria', 6 )

% 50% Acceptance criteria met
figure
hold on
% title('50% First Story Columns Acceptance Criteria') 
% ASCE 41 frags
% plot(x_points,frag_curves.io_cols_1.frag_curve_ew_50,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 IO')
% plot(x_points,frag_curves.ls_cols_1.frag_curve_ew_50,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','ASCE 41 LS')
plot(x_points,frag_curves.cp_cols_1.frag_curve_ew_50,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 CP')
% Eurcode Frags
% plot(x_points,frag_curves.euro_th_DL_cols_1.frag_curve_ew_50,'--','color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','EuroCode DL')
% plot(x_points,frag_curves.euro_th_SD_cols_1.frag_curve_ew_50,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode SD')
plot(x_points,frag_curves.euro_th_NC_cols_1.frag_curve_ew_50,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode NC')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim([0,ceil(max(sa_points))])
fn_format_and_save_plot( plot_dir, 'EuroCode - 50% Acceptance Criteria', 6 )

end

function [ ] = fn_plot_frag_curves(x_points, frag_curves, col_frag_curve, icsb_motion, plot_name, plot_dir, plot_title, x_range, matlab_colors, frag, num_gms, sa_points)    
% Import packages
import plotting_tools.fn_format_and_save_plot

hold on
sct = scatter(sa_points,frag.num_first./num_gms,[],matlab_colors(1,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Element')
sct = scatter(sa_points,frag.num_10./num_gms,[],matlab_colors(2,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Mechanism')
sct = scatter(sa_points,frag.num_25./num_gms,[],matlab_colors(3,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Mechanism')
sct = scatter(sa_points,frag.num_50./num_gms,[],matlab_colors(4,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Mechanism')
sct = scatter(sa_points,frag.num_75./num_gms,[],matlab_colors(5,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Mechanism')
sct = scatter(sa_points,frag.num_100./num_gms,[],matlab_colors(6,:),'filled');
set(get(get(sct,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
plot(x_points,frag_curves.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Mechanism')
plot(x_points,col_frag_curve,'k','lineWidth',2,'DisplayName','Collapse')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim(x_range)
fn_format_and_save_plot( plot_dir, [plot_name '-alt'], 6 )

hold on
% title(strrep(plot_title,'_',' '))
plot(x_points,frag_curves.frag_curve_ew_first,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','First Element')
plot(x_points,frag_curves.frag_curve_ew_10,'color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','10% of Mechanism')
plot(x_points,frag_curves.frag_curve_ew_25,'color',matlab_colors(3,:),'lineWidth',1.5,'DisplayName','25% of Mechanism')
plot(x_points,frag_curves.frag_curve_ew_50,'color',matlab_colors(4,:),'lineWidth',1.5,'DisplayName','50% of Mechanism')
plot(x_points,frag_curves.frag_curve_ew_75,'color',matlab_colors(5,:),'lineWidth',1.5,'DisplayName','75% of Mechanism')
plot(x_points,frag_curves.frag_curve_ew_100,'color',matlab_colors(6,:),'lineWidth',1.5,'DisplayName','100% of Mechanism')
plot(x_points,col_frag_curve,'k','lineWidth',2,'DisplayName','Collapse')
plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim(x_range)
fn_format_and_save_plot( plot_dir, plot_name, 6 )

end