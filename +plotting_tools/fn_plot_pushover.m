function [ model ] = fn_plot_pushover( read_dir, seismic_wt, base_shear, model )
% Description: Fn to plot all things pushover

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:


%% Initial Setup
% Import Packages
import plotting_tools.*

% Define Color Scheme
[matlab_colors] = fn_matlab_colormap;

% Load data
load([read_dir filesep 'node_analysis.mat'])

%% Begin Method
control_nodes = node(node.primary_story == 1,:);

% Calculate roof disp
roof_node = control_nodes(control_nodes.y == max(control_nodes.y),:);
load([read_dir filesep 'node_TH_' num2str(roof_node.id) '.mat'])
roof_drift_x = abs(nd_TH.disp_x_TH)/roof_node.y;
roof_drift_x_neg = abs(nd_TH.disp_x_neg_TH)/roof_node.y;
v_ratio_x = abs(base_shear.base_shear_x_TH)/sum(seismic_wt);
v_ratio_x_neg = abs(base_shear.base_shear_x_neg_TH)/sum(seismic_wt);
model.Vy_x = max([v_ratio_x,v_ratio_x_neg]);
if isfield(nd_TH,'disp_z_TH')
    roof_drift_z = abs(nd_TH.disp_z_TH)/roof_node.y;
    roof_drift_z_neg = abs(nd_TH.disp_z_neg_TH)/roof_node.y;
    v_ratio_z = abs(base_shear.base_shear_z_TH)/sum(seismic_wt);
    v_ratio_z_neg = abs(base_shear.base_shear_z_neg_TH)/sum(seismic_wt);
    model.Vy_z = max([v_ratio_z,v_ratio_z_neg]);
end

% Plot Roof Drift Pushover Normalized by Building Weight
hold on
if isfield(nd_TH,'disp_z_TH')
    plot(roof_drift_z,v_ratio_z,'color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','NS (Y)')
    plot(roof_drift_z_neg,v_ratio_z_neg,'--','color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','NS (-Y)')
end
plot(roof_drift_x,v_ratio_x,'color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','EW (X)')
plot(roof_drift_x_neg,v_ratio_x_neg,'--','color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','EW (-X)')
ylabel('Base Shear / Seismic Weight')
xlabel('Roof Drift')
box on
legend('location','northeast')
legend boxoff
set(gca,'FontSize',12)
fn_format_and_save_plot( [read_dir filesep 'Pushover_Plots'], 'Normalized Pushover', 0 )

end

