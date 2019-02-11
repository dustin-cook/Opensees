function [ ] = fn_plot_pushover( read_dir, direction, seismic_wt, base_shear )
% Description: Fn to plot all things pushover

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:


%% Initial Setup
% Import Packages
import plotting_tools.fn_format_and_save_plot

% Load data
load([read_dir filesep 'node_analysis.mat'])

%% Begin Method
control_nodes = node(node.primary_story == 1,:);
base_nodes = node(node.y == 0,:);

% Calculate roof disp
roof_node = control_nodes(control_nodes.y == max(control_nodes.y),:);
roof_disp = roof_node.(['disp_' direction '_TH']){1};
roof_drift = roof_disp/roof_node.y;

% Plot Roof Drift Pushover Normalized by Building Weight
v_ratio = base_shear/sum(seismic_wt);
plot(roof_drift,v_ratio)
ylabel('Base Shear / Seismic Weight')
xlabel('Roof Drift')
plot_dir = [read_dir filesep 'Pushover_Plots'];
plot_name = ['Normalized Pushover - ' direction];
fn_format_and_save_plot( plot_dir, plot_name, 2 )

% Plot story Drift Pushover
hold on
for i = 1:height(control_nodes)
    story_node = control_nodes(i,:);
    story_disp(i,:) = story_node.(['disp_' direction '_TH']){1};
    if i == 1
        rel_story_disp = story_disp;
        story_drift = rel_story_disp ./ story_node.y;
    else
        rel_story_disp = story_disp(i,:) - story_disp(i-1,:);
        story_drift = rel_story_disp ./ (control_nodes.y(i) - control_nodes.y(i-1));
    end
    plot(story_drift,base_shear/seismic_wt(i),'DisplayName',['Story - ' num2str(i)'])
end
ylabel('Total Base Shear (k)')
xlabel('Story Drift (in)')
plot_dir = [read_dir filesep 'Pushover_Plots'];
plot_name = ['Story Pushover - ' direction];
fn_format_and_save_plot( plot_dir, plot_name, 1 )
end

