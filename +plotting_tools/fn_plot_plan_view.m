function [ ] = fn_plot_plan_view( hinge, element, node, ele_side, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Begin Method
% Filter Elements
ele_2_use = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,ele_2_use.node_1),:);
hinge_2_use = hinge(ismember(hinge.element_id,ele_2_use.id) & hinge.ele_side == ele_side,:);
color_range = [0 1.2];

% %% Plot A Ratio
% fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - A ratio'], 'Max(\theta)/"a"', color_range )
% fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - A ratio OOP'], 'Max(\theta)/"a"', color_range )
% srss_value = sqrt(hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')).^2); % Update to 1.5 for linear value
% fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - A ratio SRSS'], 'Max(\theta)/"a"', color_range )
% 
% %% Plot B Ratio
% fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - B ratio'], 'Max(\theta)/"b"', color_range )
% fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - B ratio OOP'], 'Max(\theta)/"b"' , color_range)
% srss_value = sqrt(hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
% fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - B ratio SRSS'], 'Max(\theta)/"b"', color_range )
% 
% %% Plot V Ratio
% fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - V ratio'], 'Max(V)/Vn', color_range )
% fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_name, ' - V ratio OOP'], 'Max(V)/Vn', color_range )
% srss_value = sqrt(hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
% fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_name, ' - V ratio SRSS'], 'Max(V)/Vn', color_range )

% %% Plot P Ratio
% fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_expected(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - P ratio expected'], 'P_{max} / f''_{c,e}A_g', color_range )
% fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_nominal(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - P ratio nominal'], 'P_{max} / f''_{c,n}A_g', color_range )
% fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_asce41(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - P ratio asce41'], 'P_{max} / P_n', color_range )

%% Plot TAR
% story_nodes = node(node.story == 1 & node.record_accel == 1 & node.mass > 0 & node.z ~= 450,:);
% color_range = [0 2];
% fn_plot_plan_scatter( story_nodes, story_nodes.TAR_x, plot_dir, [plot_name, ' - TAR x'], 'TAR - EW Shaking', color_range )
% fn_plot_plan_scatter( story_nodes, story_nodes.TAR_z, plot_dir, [plot_name, ' - TAR z'], 'TAR - NS Shaking', color_range )
% fn_plot_plan_scatter( story_nodes, story_nodes.TAR_srss, plot_dir, [plot_name, ' - TAR srrs'], 'TAR - SRSS', color_range )

%% Plot Recorded Damage
color_range = [0 3];
fn_plot_plan_scatter( node_2_use, hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_name, ' - Recorded Damage'], 'Recorded Damage', color_range )

%% Plotter
function [ ] = fn_plot_plan_scatter( node_2_use, scatter_value, plot_dir, plot_name, plot_title, color_range )
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
caxis(color_range)
scatter(node_2_use.x,node_2_use.z,400,scatter_value,'s','filled')
colorbar
set(gca,'FontSize',15)
axis off
title(plot_title)

% Format and save plot
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end


end

