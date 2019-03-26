function [ ] = fn_plot_plan_view( hinge, element, node, ele_side, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Begin Method
% Filter Elements
ele_2_use = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,ele_2_use.node_1),:);
hinge_2_use = hinge(ismember(hinge.element_id,ele_2_use.id) & hinge.ele_side == ele_side,:);

%% Plot A Ratio
fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - A ratio'], 'Max(\theta)/"a"' )
fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - A ratio OOP'], 'Max(\theta)/"a"' )
srss_value = sqrt(hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')).^2); % Update to 1.5 for linear value
fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - A ratio SRSS'], 'Max(\theta)/"a"' )

%% Plot B Ratio
fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - B ratio'], 'Max(\theta)/"b"' )
fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - B ratio OOP'], 'Max(\theta)/"b"' )
srss_value = sqrt(hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - B ratio SRSS'], 'Max(\theta)/"b"' )

%% Plot V Ratio
fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - V ratio'], 'Max(V)/Vn' )
fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - V ratio OOP'], 'Max(V)/Vn' )
srss_value = sqrt(hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - V ratio SRSS'], 'Max(V)/Vn' )

%% Plotter
function [ ] = fn_plot_plan_scatter( node_2_use, scatter_value, plot_dir, plot_name, plot_title )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial setup
% Import Packages
import plotting_tools.*

% Define color scheme
colormap jet

% base shaded region
r_width = max(node_2_use.x)-min(node_2_use.x);
r_height = max(node_2_use.z)-min(node_2_use.z);

%% Plot
hold on
rectangle('Position',[min(node_2_use.x),min(node_2_use.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.5,0.5,0.5])
set(gca,'XTick',[],'YTick',[])
xlim([min(node_2_use.x)-r_width*0.05,max(node_2_use.x)+r_width*0.05])
ylim([min(node_2_use.z)-r_height*0.05,max(node_2_use.z)+r_height*0.05])

% Plot Column Highlight
caxis([0,1.2])
scatter(node_2_use.x,node_2_use.z,400,scatter_value,'s','filled')
colorbar
set(gca,'FontSize',15)
axis off
title(plot_title)

% Format and save plot
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end


end

