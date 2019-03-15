function [ ] = fn_plot_plan_view( hinge, element, node, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial setup
% Import Packages
import plotting_tools.*

% Define color scheme
colormap jet

%% Begin Method
% Filter Elements
ele_2_use = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,ele_2_use.node_1),:);
hinge_2_use = hinge(ismember(hinge.element_id,ele_2_use.id) & strcmp(hinge.direction,'primary') & hinge.ele_side == 1,:);

% base shaded region
r_width = max(node_2_use.x)-min(node_2_use.x);
r_height = max(node_2_use.z)-min(node_2_use.z);

%% Plot A Ratio
subplot(2,2,1)
hold on
rectangle('Position',[min(node_2_use.x),min(node_2_use.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.5,0.5,0.5])
set(gca,'XTick',[],'YTick',[])
xlim([min(node_2_use.x)-r_width*0.05,max(node_2_use.x)+r_width*0.05])
ylim([min(node_2_use.z)-r_height*0.05,max(node_2_use.z)+r_height*0.05])

% Plot Column Highlight
c = hinge_2_use.a_ratio;
caxis([0,1.2])
scatter(node_2_use.x,node_2_use.z,400,c,'s','filled')
colorbar
set(gca,'FontSize',15)
axis off
title('Max(\theta)/"a"')

%% Plot B Ratio
subplot(2,2,2)
hold on
rectangle('Position',[min(node_2_use.x),min(node_2_use.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.5,0.5,0.5])
set(gca,'XTick',[],'YTick',[])
xlim([min(node_2_use.x)-r_width*0.05,max(node_2_use.x)+r_width*0.05])
ylim([min(node_2_use.z)-r_height*0.05,max(node_2_use.z)+r_height*0.05])

% Plot Column Highlight
c = hinge_2_use.b_ratio;
caxis([0,1.2])
scatter(node_2_use.x,node_2_use.z,400,c,'s','filled')
colorbar
set(gca,'FontSize',15)
axis off
title('Max(\theta)/"b"')

%% Plot V Ratio
subplot(2,2,3)
hold on
rectangle('Position',[min(node_2_use.x),min(node_2_use.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.5,0.5,0.5])
set(gca,'XTick',[],'YTick',[])
xlim([min(node_2_use.x)-r_width*0.05,max(node_2_use.x)+r_width*0.05])
ylim([min(node_2_use.z)-r_height*0.05,max(node_2_use.z)+r_height*0.05])

% Plot Column Highlight
c = hinge_2_use.V_ratio;
caxis([0,1.2])
scatter(node_2_use.x,node_2_use.z,400,c,'s','filled')
colorbar
axis off
title('Max(V)/Vn')


% Format and save plot
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end

