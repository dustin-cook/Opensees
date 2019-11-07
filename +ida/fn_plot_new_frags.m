function [ ] = fn_plot_new_frags(analysis, model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Load data
read_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
load([read_dir filesep 'new_frag_curves.mat'])
load([read_dir filesep 'gm_data.mat'])

%% Plot Frag Curves
% color scheme 
matlab_colors = [0 0.4470 0.7410;
          0.85 0.325 0.0980;
          0.929 0.694 0.1250;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.7445 0.933;
          0.635 0.078 0.184];
matlab_colors = [matlab_colors;matlab_colors;matlab_colors]; % stack matlab colors so they repeat
      
plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'New Fragility Plots'];


%% Plot 1 - Collapse v Min Mean Max Rot Limit
hold on
max_val = 16;
x_points = linspace(max_val/100,max_val,100);

rank_val = sort(gm_data.collapse.cols_walls_1_min_b_e(~isnan(gm_data.collapse.cols_walls_1_min_b_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_min_b_e.theta),new_frag_curves.collapse.cols_walls_1_min_b_e.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','Min of Components')

rank_val = sort(gm_data.collapse.cols_walls_1_mean_b_e(~isnan(gm_data.collapse.cols_walls_1_mean_b_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_mean_b_e.theta),new_frag_curves.collapse.cols_walls_1_mean_b_e.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','Mean of Components')

rank_val = sort(gm_data.collapse.cols_walls_1_max_b_e(~isnan(gm_data.collapse.cols_walls_1_max_b_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'r','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_max_b_e.theta),new_frag_curves.collapse.cols_walls_1_max_b_e.beta);
plot(x_points,cdf,'r','lineWidth',1,'DisplayName','Max of Components')

xlabel('Deformation Demand / Deformation Capacity')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Rot Lim min-mean-max';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 2 - Collapse v Min Mean Max CP
hold on
max_val = 16;
x_points = linspace(max_val/100,max_val,100);

rank_val = sort(gm_data.collapse.cols_walls_1_min_cp(~isnan(gm_data.collapse.cols_walls_1_min_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_min_cp.theta),new_frag_curves.collapse.cols_walls_1_min_cp.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','Min of Components')

rank_val = sort(gm_data.collapse.cols_walls_1_mean_cp(~isnan(gm_data.collapse.cols_walls_1_mean_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_mean_cp.theta),new_frag_curves.collapse.cols_walls_1_mean_cp.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','Mean of Components')

rank_val = sort(gm_data.collapse.cols_walls_1_max_cp(~isnan(gm_data.collapse.cols_walls_1_max_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'r','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.cols_walls_1_max_cp.theta),new_frag_curves.collapse.cols_walls_1_max_cp.beta);
plot(x_points,cdf,'r','lineWidth',1,'DisplayName','Max of Components')

xlabel('Deformation Demand / CP Limit-State')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v CP min-mean-max';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 3 - Collapse v Percent rotation limit
max_val = 1;
hold on
rank_val = sort(gm_data.collapse.cols_walls_1_percent_b_e(~isnan(gm_data.collapse.cols_walls_1_percent_b_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1,'DisplayName','Sidesway Collapse')
xlabel('Percent of Mechanism Exceeding Rotation Limit')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Percent Rot Limit';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 4 - Collapse v Percent CP
max_val = 1;
hold on
rank_val = sort(gm_data.collapse.cols_walls_1_percent_cp(~isnan(gm_data.collapse.cols_walls_1_percent_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1,'DisplayName','Sidesway Collapse')
xlabel('Percent of Mechanism Exceeding CP')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Percent CP';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 5 - Collapse v Max Drift
max_val = 0.1;
x_points = linspace(max_val/100,max_val,100);
hold on
rank_val = sort(gm_data.collapse.max_drift(~isnan(gm_data.collapse.max_drift)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.max_drift.theta),new_frag_curves.collapse.max_drift.beta);
plot(x_points,cdf,'k','lineWidth',1)
xlabel('Max Drift')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Max Drift';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 6 - Collapse v Normalized Energy
hold on
max_val = 7;
x_points = linspace(max_val/100,max_val,100);
rank_val = sort(gm_data.collapse.norm_energy_tot(~isnan(gm_data.collapse.norm_energy_tot)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse.norm_energy_tot.theta),new_frag_curves.collapse.norm_energy_tot.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','Total Energy')
xlabel('Normalized Energy Dissapation')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Norm Energy';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 7 - Collapse v Mean Rotation Limit direction breakdown
max_val = 10;
x_points = linspace(max_val/100,max_val,100);
hold on
rank_val = sort(gm_data.collapse_x.cols_1_mean_b(~isnan(gm_data.collapse_x.cols_1_mean_b)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_x.cols_1_mean_b.theta),new_frag_curves.collapse_x.cols_1_mean_b.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','EW Frame')
rank_val = sort(gm_data.collapse_z.walls_1_mean_e(~isnan(gm_data.collapse_z.walls_1_mean_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_z.walls_1_mean_e.theta),new_frag_curves.collapse_z.walls_1_mean_e.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','NS Wall')
xlabel('Mean DCR_{Rotation Limit}')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Mean Rot Lim directions';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 8 - Collapse v Mean Rotation Limit direction breakdown
max_val = 10;
x_points = linspace(max_val/100,max_val,100);
hold on
rank_val = sort(gm_data.collapse_x.cols_1_mean_cp(~isnan(gm_data.collapse_x.cols_1_mean_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_x.cols_1_mean_cp.theta),new_frag_curves.collapse_x.cols_1_mean_cp.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','EW Frame')
rank_val = sort(gm_data.collapse_z.walls_1_mean_cp(~isnan(gm_data.collapse_z.walls_1_mean_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_z.walls_1_mean_cp.theta),new_frag_curves.collapse_z.walls_1_mean_cp.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','NS Wall')
xlabel('Mean DCR_{CP}')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Mean CP directions';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 9 - Collapse v percent Rotation Limit direction breakdown
max_val = 1;
hold on
rank_val = sort(gm_data.collapse_x.cols_1_percent_b(~isnan(gm_data.collapse_x.cols_1_percent_b)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1,'DisplayName','EW Frame')
rank_val = sort(gm_data.collapse_z.walls_1_percent_e(~isnan(gm_data.collapse_z.walls_1_percent_e)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
plot(rank_val,rank,'b','lineWidth',1,'DisplayName','NS Wall')
xlabel('Percent of Components Exceeding Rotation Limit')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Percent Rotation Limit directions';
legend('location','northwest')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 10 - Collapse v percent CP direction breakdown
max_val = 1;
hold on
rank_val = sort(gm_data.collapse_x.cols_1_percent_cp(~isnan(gm_data.collapse_x.cols_1_percent_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1,'DisplayName','EW Frame')
rank_val = sort(gm_data.collapse_z.walls_1_percent_cp(~isnan(gm_data.collapse_z.walls_1_percent_cp)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
plot(rank_val,rank,'b','lineWidth',1,'DisplayName','NS Wall')
xlabel('Percent of Components Exceeding CP')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Percent CP directions';
legend('location','northwest')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 11 - Collapse v Drift per direction
hold on
max_val = 0.1;
x_points = linspace(max_val/100,max_val,100);
rank_val = sort(gm_data.collapse_x.drift_x(~isnan(gm_data.collapse_x.drift_x)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_x.drift_x.theta),new_frag_curves.collapse_x.drift_x.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','EW Columns')
rank_val = sort(gm_data.collapse_z.drift_z(~isnan(gm_data.collapse_z.drift_z)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_z.drift_z.theta),new_frag_curves.collapse_z.drift_z.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','NS Walls')
xlabel('Max Drift')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Drfit directions';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 12 -  Collapse v Normalized Energy per direction
hold on
max_val = 10;
x_points = linspace(max_val/100,max_val,100);
rank_val = sort(gm_data.collapse_x.norm_energy_ew(~isnan(gm_data.collapse_x.norm_energy_ew)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_x.norm_energy_ew.theta),new_frag_curves.collapse_x.norm_energy_ew.beta);
plot(x_points,cdf,'k','lineWidth',1,'DisplayName','EW Columns')
rank_val = sort(gm_data.collapse_z.norm_energy_ns(~isnan(gm_data.collapse_z.norm_energy_ns)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(new_frag_curves.collapse_z.norm_energy_ns.theta),new_frag_curves.collapse_z.norm_energy_ns.beta);
plot(x_points,cdf,'b','lineWidth',1,'DisplayName','NS Walls')
xlabel('Normalized Energy Dissipation')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Norm Energy directions';
legend('location','southeast')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 13 - Collapse v gravity
max_val = 1;
hold on
rank_val = sort(gm_data.collapse.gravity_load_lost_ratio(~isnan(gm_data.collapse.gravity_load_lost_ratio)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1)
xlabel('Percent of Gravity Load Support Failure')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Gravity';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Plot 14 - Collapse v gravity direction breakdown
max_val = 1;
hold on
rank_val = sort(gm_data.collapse_x.gravity_load_lost_ratio(~isnan(gm_data.collapse_x.gravity_load_lost_ratio)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'k','filled','HandleVisibility','off');
plot(rank_val,rank,'k','lineWidth',1,'DisplayName','EW Frame')
rank_val = sort(gm_data.collapse_z.gravity_load_lost_ratio(~isnan(gm_data.collapse_z.gravity_load_lost_ratio)));
rank = (1:length(rank_val))/length(rank_val);
scatter(rank_val,rank,'b','filled','HandleVisibility','off');
plot(rank_val,rank,'b','lineWidth',1,'DisplayName','NS Wall')
xlabel('Percent of Gravity Load Support Failure')
ylabel('P[Collapse]')
xlim([0,max_val])
plot_name = 'Collapse v Gravity directions';
legend('location','northwest')
fn_format_and_save_plot( plot_dir, plot_name, 4 )


% % General Plots (ie all)
% figure
% hold on
% flds_1 = fieldnames(new_frag_curves);
% for i = 1:length(flds_1)
%     struct_name = flds_1{i};
%     flds_2 = fieldnames(new_frag_curves.(struct_name));
%     plot_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'New Fragility Plots' '/' struct_name];
%     for j = 1:length(flds_2)
%         hold on
%         var_name = flds_2{j};
%         max_val = 2*new_frag_curves.(struct_name).(var_name).theta;
%         rank_val = sort(gm_data.(struct_name).(var_name)(~isnan(gm_data.(struct_name).(var_name))));
%         rank = (1:length(rank_val))/length(rank_val);
%         scatter(rank_val,rank,'k','filled','HandleVisibility','off');
%         x_points = linspace(max_val/100,max_val,100);
%         cdf = logncdf(x_points,log(new_frag_curves.(struct_name).(var_name).theta),new_frag_curves.(struct_name).(var_name).beta);
%         plot(x_points,cdf,'k','lineWidth',1)
%         xlabel(strrep(var_name,'_',' '))
%         ylabel(strrep(struct_name,'_',' '))
%         xlim([0,max_val])
%         plot_name = [strrep(struct_name,'_',' ') ' v ' strrep(var_name,'_',' ')];
%         fn_format_and_save_plot( plot_dir, plot_name, 4 )
%     end
% end
 