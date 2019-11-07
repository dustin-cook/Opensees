function [ ] = fn_plot_plan_scatter( node_2_use, scatter_value, plot_dir, plot_name, discrete, plot_txt, x_lab, recorded_damage, c_range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial setup
% Import Packages
import plotting_tools.*

% base shaded region
r_width = max(node_2_use.x)-min(node_2_use.x);
r_height = max(node_2_use.z)-min(node_2_use.z);
x_center = (max(node_2_use.x)-min(node_2_use.x))/2 + min(node_2_use.x);
z_center = (max(node_2_use.z)-min(node_2_use.z))/2 + min(node_2_use.z);

%% Plot
hold on
rectangle('Position',[min(node_2_use.x),min(node_2_use.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.15,0.15,0.15],'LineWidth',1.5)
set(gca,'XTick',[],'YTick',[])
plot_width = max([r_height,r_width])*1.05;
xlim([x_center-plot_width/2,x_center+plot_width/2])
ylim([z_center-plot_width/2,z_center+plot_width/2])

% Plot Column Highlight
% caxis(color_range)
% Define color scheme
if discrete == 3
    cmap = [145,207,96;
            252,141,89]/255;
    colormap(cmap);
    caxis([0 1]) 
elseif discrete
    cmap = [200 200 200;
            238 235 112;
            252 160 98;
            241 95 95]/255;
    colormap(cmap);
    caxis([1 4])
else
    div_colors = [0,60,48;
              1,102,94;
              53,151,143;
              128,205,193;
              199,234,229;
              246,232,195;
              223,194,125;
              191,129,45;
              140,81,10;
              84,48,5]/255;
          
    colormap(div_colors)
    caxis([0 c_range])
end
scatter(node_2_use.x,node_2_use.z,150,scatter_value,'s','filled','LineWidth',1.25,'MarkerEdgeColor',[0.15 0.15 0.15])
% text(x_center,z_center,plot_txt,...
%     'FontWeight','bold',...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle')
axis off

if discrete == 0
    c = colorbar('south');
    c.AxisLocation = 'out';
    c.Label.String = x_lab;
%     scatter(node_2_use.x(recorded_damage >= 2),node_2_use.z(recorded_damage >= 2),350,'k','s','LineWidth',1.15)
    set(gca,'FontSize',14)
elseif discrete == 2
    c = colorbar('south');
    c.AxisLocation = 'out';
    c.Ticks = [0.4 1.1 1.9 2.6];
    c.TickLabels = {'IO', 'LS', 'CP', 'Fail All'};
    c.Label.String = x_lab;
%     scatter(node_2_use.x(recorded_damage >= 2),node_2_use.z(recorded_damage >= 2),350,'k','s','LineWidth',1.15)
    set(gca,'FontSize',14)
elseif discrete == 3
%     scatter(node_2_use.x(recorded_damage >= 2),node_2_use.z(recorded_damage >= 2),350,'k','s','LineWidth',1.15)
end
% Format and save plot
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end