function [ model ] = fn_plot_pushover( read_dir, story, base_shear, model )
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
v_ratio_x = abs(base_shear.base_shear_x_TH)/sum(story.seismic_wt);
v_ratio_x_neg = abs(base_shear.base_shear_x_neg_TH)/sum(story.seismic_wt);
model.Vy_x = max([v_ratio_x,v_ratio_x_neg]);
if isfield(nd_TH,'disp_z_TH')
    roof_drift_z = abs(nd_TH.disp_z_TH)/roof_node.y;
    roof_drift_z_neg = abs(nd_TH.disp_z_neg_TH)/roof_node.y;
    v_ratio_z = abs(base_shear.base_shear_z_TH)/sum(story.seismic_wt);
    v_ratio_z_neg = abs(base_shear.base_shear_z_neg_TH)/sum(story.seismic_wt);
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

% Plot Roof Drift Pushover vs Baseshear (kips)
hold on
if isfield(nd_TH,'disp_z_TH')
    plot(roof_drift_z,abs(base_shear.base_shear_z_TH)/1000,'color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','NS (Y)')
    plot(roof_drift_z_neg,abs(base_shear.base_shear_z_neg_TH)/1000,'--','color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','NS (-Y)')
end
plot(roof_drift_x,abs(base_shear.base_shear_x_TH)/1000,'color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','EW (X)')
plot(roof_drift_x_neg,abs(base_shear.base_shear_x_neg_TH)/1000,'--','color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','EW (-X)')
ylabel('Base Shear (kips)')
xlabel('Roof Drift')
box on
legend('location','northeast')
legend boxoff
set(gca,'FontSize',12)
fn_format_and_save_plot( [read_dir filesep 'Pushover_Plots'], 'BaseShear Pushover', 0 )

% Calculate Peak drfits from pushover
for n = 1:height(control_nodes)
    story_control_node = control_nodes(n,:);
    load([read_dir filesep 'node_TH_' num2str(story_control_node.id) '.mat'])
    story.story_disp_x(story.id == story_control_node.story) = max(abs([nd_TH.disp_x_TH,nd_TH.disp_x_neg_TH]));
    if isfield(nd_TH,'disp_z_TH')
        story.story_disp_z(story.id == story_control_node.story) = max(abs([nd_TH.disp_z_TH,nd_TH.disp_z_neg_TH]));
    end
end

story.rel_disp_x = story.story_disp_x - [0; story.story_disp_x(1:(end-1))];
story.drift_x = story.rel_disp_x ./ story.story_ht;

% Plot pushover drift profile
story_vector = sort([story.id; story.id(2:end); story.id(end)+1]);
tmp_array(:,1) = [story.id; story.id];
tmp_array(:,2) = [story.drift_x; story.drift_x];
% story_vector_tab = struct2table(story_vector);
tmp_array = sortrows(tmp_array);
drift_vector = tmp_array(:,2);

plot(drift_vector,story_vector,'linewidth',1.5)
ylabel('Floor')
xlabel('Peak Drift')
box on
legend boxoff
set(gca,'FontSize',12)
fn_format_and_save_plot( [read_dir filesep 'Pushover_Plots'], 'Drift Profile', 0 )

end

