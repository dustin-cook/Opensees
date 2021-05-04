function [ ] = fn_plot_ida(analysis, model, gm_set_table, ida_results, SSF_ew, main_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
if analysis.run_z_motion 
    params = {'b','e','cp'};
else
    params = {'b','cp'};
end
% params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
frag_probs = [10 25 50 75 100];

% Define Directiries
read_dir = [main_dir '/' 'IDA' '/' 'Fragility Data'];
plot_dir = [main_dir '/' 'IDA' '/' 'IDA Plots'];
if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end

% Load data
ida_table = readtable([read_dir filesep 'ida_table.csv']);
gm_table = readtable([read_dir filesep 'gm_table.csv']);
load([read_dir filesep 'frag_curves.mat'])
% frag_curves.global = readtable([read_dir filesep 'global_limit_states.csv']);
% frag_curves.local  = readtable([read_dir filesep 'local_limit_states.csv']);

x_points = 0.1:0.01:3;

%% Plot IDA curves
% fprintf('Saving IDA Summary Data and Figures to Directory: %s \n',plot_dir)
% if analysis.run_z_motion
%     drs = {'x','z'};    
% else
%     drs = {'x'};
% end
% for d = 1:length(drs)
%     % Plot IDA curves
%     hold on
%     for gms = 1:height(gm_set_table)
%         ida_plt = plot(ida_table.(['drift_' drs{d}])(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),ida_table.(['sa_' drs{d}])(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),'-o','color',[0.75 0.75 0.75],'HandleVisibility','off');
%     end
%     % plot(frag.mean_idr_ew,frag.p695_sa_ew,'b','lineWidth',1.5,'DisplayName','Mean Drift')
%     % plot(frag.mean_idr_ew,frag.p_15_sa_ew,'--b','lineWidth',1.5,'DisplayName','15th Percentile')
%     % plot(frag.mean_idr_ew,frag.p_85_sa_ew,'--b','lineWidth',1.5,'DisplayName','85th Percentile')
%     plot([0,0.1],[ida_results.spectra(d),ida_results.spectra(d)],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
%     xlim([0 0.1])
%     xlabel('Max Drift')
%     ylabel(['Sa(T_1=' num2str(ida_results.period(d)) 's) (g)'])
%     fn_format_and_save_plot( plot_dir, ['IDA Plot ' drs{d} ' Frame Direction'], 3 )
% end

%% Component IDA Curves (only for x dir)
% vals2plot = {'gravity_dcr', 'mean_b', 'max_b', 'max_cp', 'lat_cap_ratio_both'};
% for v = 1:length(vals2plot)
%     % replace drift
%     hold on
%     for gms = 1:height(gm_set_table)
%         ida_plt = plot(ida_table.(vals2plot{v})(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),ida_table.('sa_x')(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),'-o','color',[0.75 0.75 0.75],'HandleVisibility','off');
%     end
%     xlabel(vals2plot{v})
%     ylabel(['Sa(T_1=' num2str(ida_results.period(d)) 's) (g)'])
%     fn_format_and_save_plot( plot_dir, ['IDA Plot - ' vals2plot{v} ' - 1'], 3 )
%     
% %     % replace sa
% %     hold on
% %     for gms = 1:height(gm_set_table)
% %         ida_plt = plot(ida_table.('drift_x')(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),ida_table.(vals2plot{v})(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),'-o','color',[0.75 0.75 0.75],'HandleVisibility','off');
% %     end
% %     xlabel('Max Drift')
% %     ylabel(vals2plot{v})
% %     fn_format_and_save_plot( plot_dir, ['IDA Plot - ' vals2plot{v} ' - 2'], 3 )
% end

%% Calculate Post Fragility Curve P695 factors
% if analysis.run_z_motion
%     factor_3D = 1.11;
% else
%     factor_3D = 1.0;
% end
% SSF = SSF_ew; % Assume EW SSF
% 
% % Set Ida summary table 
% ida_summary_table.direction = 'EW';
% ida_summary_table.period = ida_results.period(1);
% ida_summary_table.spectra = ida_results.spectra(1);
% ida_summary_table.mce = ida_results.mce(1);
% 
% % Collapse Median Sa
% ida_summary_table.sa_med_col = frag_curves.collapse.theta;
% ida_summary_table.sa_med_UR = frag_curves.UR.theta;
% ida_summary_table.sa_med_CP = frag_curves.b.theta(1);
% 
% % Collapse Margin Ratio
% ida_summary_table.cmr = ida_summary_table.sa_med_col / ida_summary_table.mce;
% ida_summary_table.cmr_UR = ida_summary_table.sa_med_UR / ida_summary_table.mce;
% ida_summary_table.cmr_CP = ida_summary_table.sa_med_CP / ida_summary_table.mce;
% 
% % Adjust for SSF and 3D
% ida_summary_table.acmr = factor_3D*SSF*ida_summary_table.cmr;
% ida_summary_table.acmr_UR = factor_3D*SSF*ida_summary_table.cmr_UR;
% ida_summary_table.acmr_CP = factor_3D*SSF*ida_summary_table.cmr_CP;
% median_adjustment = ida_summary_table.sa_med_col*(factor_3D*SSF - 1);
% 
% ida_summary_table.p_col_mce = logncdf(ida_summary_table.mce,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% ida_summary_table.p_col_mce_adjust = logncdf(ida_summary_table.mce,log(frag_curves.collapse.theta + median_adjustment),0.6);
% ida_summary_table.p_UR_mce = logncdf(ida_summary_table.mce,log(frag_curves.UR.theta),frag_curves.UR.beta);
% ida_summary_table.p_UR_mce_adjusted = logncdf(ida_summary_table.mce,log(frag_curves.UR.theta + median_adjustment),0.6);
% ida_summary_table.p_CP_mce = logncdf(ida_summary_table.mce,log(frag_curves.b.theta(1)),frag_curves.b.beta(1));
% ida_summary_table.p_C_ICSM = logncdf(ida_summary_table.spectra,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% 
% % Save IDA results as table
% ida_results_table = struct2table(ida_summary_table);
% writetable(ida_results_table,[plot_dir filesep 'ida_results.csv'])

%% Plot Frag Curves
% % color scheme 
% matlab_colors = [0 0.4470 0.7410;
%           0.85 0.325 0.0980;
%           0.929 0.694 0.1250;
%           0.494 0.184 0.556;
%           0.466 0.674 0.188;
%           0.301 0.7445 0.933;
%           0.635 0.078 0.184];
% matlab_colors = [matlab_colors;matlab_colors;matlab_colors]; % stack matlab colors so they repeat
%       
% % Sa points to plot fragility curves
% x_points = 0.01:0.01:4;
%       
% % Total Collapse with discrete scatter
% if analysis.run_z_motion
%     opts = {'', '_x', '_z'};    
% else
%     opts = {''};
% end
% for c = 1:length(opts)
%     hold on
%     
%     % Gravity
%     rank_sa = sort(gm_table.sa_collapse_grav(~isnan(gm_table.sa_collapse_grav)));
%     rank_val = 1:length(rank_sa);
%     sct = scatter(rank_sa,rank_val./length(rank_sa),'k','filled','HandleVisibility','off');
%     cdf = logncdf(x_points,log(frag_curves.collapse_grav.theta),frag_curves.collapse_grav.beta);
%     plot(x_points,cdf,'k','DisplayName','Gravity Collapse')
%     
%     % Sidesway
%     rank_sa = sort(gm_table.(['sa_collapse' opts{c}])(~isnan(gm_table.(['sa_collapse' opts{c}]))));
%     rank_val = 1:length(rank_sa);
%     sct = scatter(rank_sa,rank_val./length(rank_sa),'m','filled','HandleVisibility','off');
%     cdf = logncdf(x_points,log(frag_curves.(['collapse' opts{c}]).theta),frag_curves.(['collapse' opts{c}]).beta);
%     plot(x_points,cdf,'m','DisplayName','Sidesway Collapse')
%     % Drift
% %     rank_sa = sort(gm_table.sa_collapse_drift_6(~isnan(gm_table.sa_collapse_drift_6)));
% %     rank_val = 1:length(rank_sa);
% %     sct = scatter(rank_sa,rank_val./length(rank_sa),'r','filled','HandleVisibility','off');
% %     cdf = logncdf(x_points,log(frag_curves.collapse_drift_6.theta),frag_curves.collapse_drift_6.beta);
% %     plot(x_points,cdf,'r','DisplayName','Collapse at 6% Drift')
% 
% %     % All
% %     rank_sa = sort(gm_table.(['sa_collapse_2' opts{c}])(~isnan(gm_table.(['sa_collapse_2' opts{c}]))));
% %     rank_val = 1:length(rank_sa);
% %     sct = scatter(rank_sa,rank_val./length(rank_sa),'k','filled','HandleVisibility','off');
% %     cdf = logncdf(x_points,log(frag_curves.(['collapse_2' opts{c}]).theta),frag_curves.(['collapse_2' opts{c}]).beta);
% %     plot(x_points,cdf,'k','DisplayName','Collapse Fragility All')
% 
%     xlabel('Sa(T_{1-EW}) (g)')
%     ylabel('P[Collapse]')
%     xlim([0,1.5])
%     fn_format_and_save_plot( plot_dir, ['Collapse Fragility' opts{c}], 6 )
% end

% % Total Collapse with 695 adjustments
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'b','lineWidth',1,'DisplayName','Collapse Fragility')
% cdf = logncdf(x_points,log(frag_curves.collapse.theta + median_adjustment),frag_curves.collapse.beta);
% plot(x_points,cdf,'--b','lineWidth',1,'DisplayName','Adjustment')
% cdf = logncdf(x_points,log(frag_curves.collapse.theta + median_adjustment),0.6);
% plot(x_points,cdf,':b','lineWidth',1.25,'DisplayName','With Uncertainty')
% plot([ida_results.mce(1),ida_results.mce(1)],[0,1],'--k','lineWidth',0.5,'DisplayName','MCE_R')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Collapse]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Collapse Fragility with 695 Uncertainty', 6 )

% % Total Collapse with breakdown in each direction
% figure
% hold on
% rank_sa = sort(gm_table.sa_collapse(~isnan(gm_table.sa_collapse)));
% rank_val = (1:length(rank_sa))/length(rank_sa);
% scatter(rank_sa,rank_val,'k','filled','HandleVisibility','off');
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'k','DisplayName','Total Collapse Fragility')
% rank_sa = sort(gm_table.sa_collapse_x(~isnan(gm_table.sa_collapse_x)));
% rank_val = (1:length(rank_sa))/length(rank_sa);
% scatter(rank_sa,rank_val,'b','filled','HandleVisibility','off');
% cdf = logncdf(x_points,log(frag_curves.ew_collapse.theta),frag_curves.ew_collapse.beta);
% plot(x_points,cdf,'b','DisplayName','EW Frame Fragility')
% rank_sa = sort(gm_table.sa_collapse_z(~isnan(gm_table.sa_collapse_z)));
% rank_val = (1:length(rank_sa))/length(rank_sa);
% scatter(rank_sa,rank_val,'r','filled','HandleVisibility','off');
% cdf = logncdf(x_points,log(frag_curves.ns_collapse.theta),frag_curves.ns_collapse.beta);
% plot(x_points,cdf,'r','DisplayName','NS Wall Fragility')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Collapse]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Collapse Fragility Both Directions', 6 )
% 
% % ASCE 41 Unacceptable Response v Total Collapse with 695 adjustments
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'b','DisplayName','Collapse Fragility')
% cdf = logncdf(x_points,log(frag_curves.UR.theta),frag_curves.UR.beta);
% plot(x_points,cdf,'-.r','DisplayName','Unacceptable Response')
% plot([ida_results.mce(1),ida_results.mce(1)],[0,1],'--k','lineWidth',0.5,'DisplayName','MCE_R')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Unacceptable Response Fragility', 6 )

% % Plot Each acceptance Criteria
% for p = 1:length(params)
%     plot_name = [params{p} ' Fragilities'];
%     fn_plot_frag_curves(x_points, frag_curves.(params{p}), frag_curves.collapse_2, frag_curves.collapse, ida_results.spectra(1), plot_name, plot_dir, [0,1.5], params{p}  )
%     if analysis.run_z_motion
%         plot_name = [params{p} ' Fragilities_x'];
%         fn_plot_frag_curves(x_points, frag_curves.([params{p} '_x']), frag_curves.collapse_2_x, frag_curves.collapse_x, ida_results.spectra(1), plot_name, plot_dir, [0,1.5], params{p}  )
%         plot_name = [params{p} ' Fragilities_z'];
%         fn_plot_frag_curves(x_points, frag_curves.([params{p} '_z']), frag_curves.collapse_2_z, frag_curves.collapse_z, ida_results.spectra(2), plot_name, plot_dir, [0,1.5], params{p}  )
%     end
% end

% First Acceptance criteria met
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.cols_1.cp.theta(1)),frag_curves.cols_1.cp.beta(1));
% plot(x_points,cdf,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 CP')
% cdf = logncdf(x_points,log(frag_curves.cols_1.euro_th_NC.theta(1)),frag_curves.cols_1.euro_th_NC.beta(1));
% plot(x_points,cdf,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode NC')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'EuroCode - First Acceptance Criteria', 6 )
% 
% % 50% Acceptance criteria met
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.cols_1.cp.theta(4)),frag_curves.cols_1.cp.beta(4));
% plot(x_points,cdf,'color',matlab_colors(1,:),'lineWidth',1.5,'DisplayName','ASCE 41 CP')
% cdf = logncdf(x_points,log(frag_curves.cols_1.euro_th_NC.theta(4)),frag_curves.cols_1.euro_th_NC.beta(4));
% plot(x_points,cdf,'--','color',matlab_colors(2,:),'lineWidth',1.5,'DisplayName','EuroCode NC')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'EuroCode - 50% Acceptance Criteria', 6 )

% % Plot For Gravity Load Lost
% figure
% hold on
% for f = 1:length(frag_probs) 
%     cdf = logncdf(x_points,log(frag_curves.gravity.(['percent_lost_' num2str(frag_probs(f))]).theta),frag_curves.gravity.(['percent_lost_' num2str(frag_probs(f))]).beta);
%     plot(x_points,cdf,'color',matlab_colors(f+1,:),'lineWidth',1.5,'DisplayName',[num2str(frag_probs(f)) '% Gravity Capacity Lost'])
% end
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'-k','lineWidth',2,'DisplayName','Collapse Fragility')
% cdf = logncdf(x_points,log(frag_curves.UR.theta),frag_curves.UR.beta);
% plot(x_points,cdf,'-.r','lineWidth',2,'DisplayName','Unacceptable Response')
% plot([ida_results.spectra(1) ida_results.spectra(1)],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Gravity Load Lost Fragilities', 6 )

% Plot Local v global fragilities
col_met = {'collapse_drift','collapse_grav','collapse_grav_or_lat'};
col_met_lin = {'--',':','-'};
local_lim_states = unique(frag_curves.local.limit_state_param);
% local_lim_state_names = strrep(strrep(local_lim_states,'*',''),'/','_');
for p = 1:length(local_lim_states)
    this_dir = [main_dir filesep 'IDA' filesep 'atc138-5' filesep local_lim_states{p}];
    mkdir(this_dir)
    
    % Initial Setup
    local_table = frag_curves.local(strcmp(frag_curves.local.limit_state_param,local_lim_states{p}),:);
    global_table = frag_curves.global(ismember(frag_curves.global.limit_state,col_met),:);
    CIR.limit_state_degree = local_table.limit_state_degree;
    RI.limit_state_degree = local_table.limit_state_degree;
    CIR_2 = {};
    num_params = height(local_table);
    if num_params > 1
        colormap(parula)
        caxis([0 num_params])
        cmap = parula(num_params);
    else
        cmap = 'b';
    end
    
    % Continuous Plots
    hold on
    % Local Values
    cdf_loc = [];
    for loc = 1:num_params
        cdf_params = local_table(loc,:);
        cdf_loc(loc,:) = logncdf(x_points,log(cdf_params.median),cdf_params.beta);
        plot(x_points,cdf_loc(loc,:),'color',cmap(loc,:),'lineWidth',1.5,'DisplayName',strrep(cdf_params.limit_state_degree{1},'_',' '))
    end

    % Global Values
    for glo = 1:height(global_table)
        cdf_params = global_table(glo,:);
        cdf_glob = logncdf(x_points,log(cdf_params.median),cdf_params.beta);
        plot(x_points,cdf_glob,[col_met_lin{glo} 'k'],'lineWidth',2,'DisplayName',strrep(cdf_params.limit_state{1},'_',' '))
        
        % Create CIR Table vals
        CIR.(cdf_params.limit_state{1}) = cdf_params.median ./ local_table.median;
        
        % Create RI Table vals
        RI.(cdf_params.limit_state{1}) = local_table.(['RI_' cdf_params.limit_state{1}]);
        
        % Save CIR2 values (P[COL|LS])
        CIR_2{glo} = cdf_glob ./ cdf_loc;
    end

    title(strrep(local_lim_states{p},'_',' '))
    xlabel('Sa(T_{1-EW}) (g)')
    ylabel('P[Exceedance]')
    xlim([0,2])
    fn_format_and_save_plot( this_dir, ['Fragility Plot - ' local_lim_states{p}], 1 )
    
    % Discrete Plots
    hold on
    % Local Values
    cdf_loc = [];
    for loc = 1:num_params
        sa_data = sort(gm_table.(['sa_' local_table.limit_state_lab{loc}]));
        rank = (1:length(sa_data))/length(sa_data);
        plot(sa_data,rank,'color',cmap(loc,:),'marker','o','lineWidth',1,'DisplayName',strrep(local_table.limit_state_degree{loc},'_',' '))
    end

    % Global Values
    for glo = 1:height(global_table)
        sa_data = sort(gm_table.(['sa_' global_table.limit_state{glo}]));
        rank = (1:length(sa_data))/length(sa_data);
        plot(sa_data,rank,[col_met_lin{glo} 'k'],'lineWidth',2,'DisplayName',strrep(global_table.limit_state{glo},'_',' '))
    end

    title(strrep(local_lim_states{p},'_',' '))
    xlabel('Sa(T_{1-EW}) (g)')
    ylabel('P[Exceedance]')
    xlim([0,2])
    fn_format_and_save_plot( this_dir, ['Fragility Plot - discrete - ' local_lim_states{p}], 1 )
    
   % Create CIR Table
   writetable(struct2table(CIR),[this_dir filesep 'cir_table.csv']);
   
   % Create CIR Table
   writetable(struct2table(RI),[this_dir filesep 'RI_table.csv']);
   
   % Create CIR 2a plot
   for glo = 1:height(global_table)
       hold on 
       for loc = 1:num_params
           plot(x_points, CIR_2{glo}(loc,:),'color',cmap(loc,:),'lineWidth',1.5,'DisplayName',strrep(local_table.limit_state_degree{loc},'_',' '))
       end
       col_title = strrep(strrep(global_table.limit_state{glo},'collapse_','collapse: '),'_',' ');
       title([strrep(local_lim_states{p},'_',' ') ' - ' col_title])
       xlabel('Sa(T_{1-EW}) (g)')
       ylabel('P[Global|Local]')
       ylim([0,1])
       xlim([0,2])
       fn_format_and_save_plot( this_dir, ['CIR 2a Plot - ' local_lim_states{p} ' - ' global_table.limit_state{glo}], 1 )
   end
   
   % Create CIR2b plot
  for glo = 1:height(global_table)
       hold on 
       for loc = 1:num_params
           cdf_params = local_table(loc,:);
           plot(cdf_params.sa_vals{1}, cdf_params.(['CIR2b_' global_table.limit_state{glo}]){1},'color',cmap(loc,:),'lineWidth',1.5,'DisplayName',strrep(local_table.limit_state_degree{loc},'_',' '))
       end
       col_title = strrep(strrep(global_table.limit_state{glo},'collapse_','collapse: '),'_',' ');
       title([strrep(local_lim_states{p},'_',' ') ' - ' col_title])
       xlabel('Sa(T_{1-EW}) (g)')
       ylabel('P[Global|Local]')
       ylim([0,1])
       xlim([0,2])
       fn_format_and_save_plot( this_dir, ['CIR 2b Plot - ' local_lim_states{p} ' - ' global_table.limit_state{glo}], 1 )
  end
   
    % CIR3 Plot
    for glo = 1:height(global_table)
       data = [];
       for loc = 1:num_params
           data(:,loc) = gm_table.(['CIR3_' local_table.limit_state_lab{loc} '_' global_table.limit_state{glo}]);
       end
       boxplot(data,'Labels',local_table.limit_state_degree)
       col_title = strrep(strrep(global_table.limit_state{glo},'collapse_','collapse: '),'_',' ');
       title([strrep(local_lim_states{p},'_',' ') ' - ' col_title])
       xlabel(['Limit State: ' strrep(local_lim_states{p},'_',' ')])
       ylabel('Collapse Margin')
       fn_format_and_save_plot( this_dir, ['CIR 3 Plot - ' local_lim_states{p} ' - ' global_table.limit_state{glo}], 0 )
    end
   
    % CIR2c Plot
    for glo = 1:height(global_table)
       data = [];
       for loc = 1:num_params
           data(:,loc) = gm_table.(['CIR2c_' local_table.limit_state_lab{loc} '_' global_table.limit_state{glo}]);
       end
       boxplot(data,'Labels',local_table.limit_state_degree)
       col_title = strrep(strrep(global_table.limit_state{glo},'collapse_','collapse: '),'_',' ');
       title([strrep(local_lim_states{p},'_',' ') ' - ' col_title])
       xlabel(['Limit State: ' strrep(local_lim_states{p},'_',' ')])
       ylabel('P[ Collapse | Limit State ]')
       fn_format_and_save_plot( this_dir, ['CIR 2c Plot - ' local_lim_states{p} ' - ' global_table.limit_state{glo}], 0 )
   end
end

% Plot For Drift Fragilities
% colormap(parula)
% caxis([0 10])
% cmap = parula(10);
% % for c = 1:length(opts)
%     hold on
% %     for d = 1:10
% 
%         % Local Values
%         cdf_params = frag_curves.local(1,:);
%         cdf = logncdf(x_points,log(cdf_params.median),cdf_params.beta);
%         plot(x_points,cdf,'color',cmap(1,:),'lineWidth',1.5,'HandleVisibility','off')
% %     end
% 
%     % Global Values
%     cdf_params = frag_curves.global(1,:);
%     cdf = logncdf(x_points,log(cdf_params.median),cdf_params.beta);
%     plot(x_points,cdf,'--k','lineWidth',2,'DisplayName',cdf_params.limit_state{1})
%     
%     xlabel('Sa(T_{1-EW}) (g)')
%     ylabel('P[Exceedance]')
%     xlim([0,2])
% %     h = colorbar;
%     ylabel(h, 'Peak Drift')
%     for i = 1:length(h.TickLabels)
%         h.TickLabels{i} = [num2str(str2double(h.TickLabels{i})*20) '%'];
%     end
%     fn_format_and_save_plot( plot_dir, ['Max Drift Fragilities' opts{c}], 6 )

%     % Plot For grav dcr
%     hold on
%     for d = 1:10
%         cdf = logncdf(x_points,log(frag_curves.gravity.(['dcr_' num2str(d*10)]).theta),frag_curves.gravity.(['dcr_' num2str(d*10)]).beta);
%         plot(x_points,cdf,'color',matlab_colors(d,:),'lineWidth',1.5,'DisplayName',['Grav DCR > ' num2str(d)/10])
%     end
%     cdf = logncdf(x_points,log(frag_curves.(['collapse'  opts{c}]).theta),frag_curves.(['collapse'  opts{c}]).beta);
%     plot(x_points,cdf,'-k','lineWidth',2,'DisplayName','Collapse Fragility')
%     xlabel('Sa(T_{1-EW}) (g)')
%     ylabel('P[Exceedance]')
%     xlim([0,2])
%     fn_format_and_save_plot( plot_dir, ['Grav DCR Fragilities' opts{c}], 6 )
% end

% % Plot For lateral cap
% hold on
% for d = 1:9
%     cdf = logncdf(x_points,log(frag_curves.lateral.(['cap_both_' num2str(d*10)]).theta),frag_curves.lateral.(['cap_both_' num2str(d*10)]).beta);
%     plot(x_points,cdf,'color',matlab_colors(d,:),'lineWidth',1.5,'DisplayName',[num2str(d)/10 '% Lat Capacity Remains'])
% end
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'-k','lineWidth',2,'DisplayName','Collapse Fragility')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Lateral Capacity Fragilities', 6 )
% 
% 
% % Plot For Adjacent Components
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.adjacent_comp.any.theta),frag_curves.adjacent_comp.any.beta);
% plot(x_points,cdf,'color',matlab_colors(1,:),'DisplayName','Any Adjacent Component')
% cdf = logncdf(x_points,log(frag_curves.adjacent_comp.any_frame.theta),frag_curves.adjacent_comp.any_frame.beta);
% plot(x_points,cdf,'color',matlab_colors(2,:),'DisplayName','Any Adjacent Frame Component')
% cdf = logncdf(x_points,log(frag_curves.adjacent_comp.all.theta),frag_curves.adjacent_comp.all.beta);
% plot(x_points,cdf,'color',matlab_colors(3,:),'DisplayName','All Adjacent Components')
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'-k','DisplayName','Collapse Fragility')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Adjacent Component Fragilities', 6 )
% 
% % Plot For Percent of Collapse Energy
% figure
% hold on
% for f = 1:length(frag_probs) 
%     cdf = logncdf(x_points,log(frag_curves.energy.(['percent_collapse_' num2str(frag_probs(f))]).theta),frag_curves.energy.(['percent_collapse_' num2str(frag_probs(f))]).beta);
%     plot(x_points,cdf,'color',matlab_colors(f,:),'DisplayName',[num2str(frag_probs(f)) '% of Collapse Energy'])
% end
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'-k','DisplayName','Collapse Fragility')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Collapse Energy Fragilities', 6 )
% 
% % Plot Multiple trigger types together
% % Plot For Adjacent Components
% figure
% hold on
% cdf = logncdf(x_points,log(frag_curves.cols_walls_1.b_e.theta(1)),frag_curves.cols_walls_1.b_e.beta(1));
% plot(x_points,cdf,'color',matlab_colors(1,:),'DisplayName','ASCE 41 CP')
% cdf = logncdf(x_points,log(frag_curves.adjacent_comp.any.theta),frag_curves.adjacent_comp.any.beta);
% plot(x_points,cdf,'color',matlab_colors(2,:),'DisplayName','Any Adjacent Component')
% cdf = logncdf(x_points,log(frag_curves.gravity.percent_lost_25.theta),frag_curves.gravity.percent_lost_25.beta);
% plot(x_points,cdf,'color',matlab_colors(3,:),'DisplayName','25% Load Capacity Lost')
% cdf = logncdf(x_points,log(frag_curves.drift.idr_2.theta),frag_curves.drift.idr_2.beta);
% plot(x_points,cdf,'color',matlab_colors(4,:),'DisplayName','2% Max Interstory Drift')
% cdf = logncdf(x_points,log(frag_curves.collapse.theta),frag_curves.collapse.beta);
% plot(x_points,cdf,'-k','DisplayName','Collapse Fragility')
% xlabel('Sa(T_{1-EW}) (g)')
% ylabel('P[Exceedance]')
% xlim([0,2])
% fn_format_and_save_plot( plot_dir, 'Multi Fragility Compare', 6 )

end

function [ ] = fn_plot_frag_curves(x_points, frag_curves, col_frag_curve, col_no_grav, icsb_motion, plot_name, plot_dir, x_range, param_name)    
% Import packages
import plotting_tools.*

% colormap(parula)
colormap(parula)
cmap = parula(100);

hold on
cdf = logncdf(x_points,log(frag_curves.theta(1)),frag_curves.beta(1));
plot(x_points,cdf,'color',cmap(round(100*frag_curves.prct_mech(1)),:),'lineWidth',1.5,'HandleVisibility','off')
for i = 2:height(frag_curves)
    cdf = logncdf(x_points,log(frag_curves.theta(i)),frag_curves.beta(i));
%     plot(x_points,cdf,'color',cmap(i,:),'lineWidth',1.5,'DisplayName',[num2str(100*frag_curves.prct_mech(i)) '% of Mechanism'])
    plot(x_points,cdf,'color',cmap(round(100*frag_curves.prct_mech(i)),:),'lineWidth',1.5,'HandleVisibility','off')
end
cdf = logncdf(x_points,log(col_frag_curve.theta),col_frag_curve.beta);
plot(x_points,cdf,'--k','lineWidth',2,'DisplayName','Gravity Collapse')
cdf = logncdf(x_points,log(col_no_grav.theta),col_no_grav.beta);
plot(x_points,cdf,'--m','lineWidth',2,'DisplayName','Sidesway Collapse')
% plot([icsb_motion icsb_motion],[0 1],'--k','lineWidth',1.5,'DisplayName','ICSB Motion')
xlabel('Sa(T_{1-EW}) (g)')
ylabel('P[Exceedance]')
xlim(x_range)
h = colorbar;
ylabel(h, ['Fraction of Components Exceeding ' upper(param_name)])
fn_format_and_save_plot( plot_dir, plot_name, 6 )

end