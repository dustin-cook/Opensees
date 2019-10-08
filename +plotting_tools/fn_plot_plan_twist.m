function [ ] = fn_plot_plan_twist( node, read_dir_opensees, plot_dir, plot_name, record_edp )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial setup
% Import Packages
import plotting_tools.*

% Hard code to ICSB for now
x_min = 71; 
x_max = 1571; 

% scale plot amplification
x_plt_scale = 100;
z_plt_scale = 100;

% Set corner cordinates
node_SW = node(node.x == x_min & node.z == min(node.z) & node.y == max(node.y),:);
node_NW = node(node.x == x_min & node.z == max(node.z) & node.y == max(node.y),:);
node_NE = node(node.x == x_max & node.z == max(node.z) & node.y == max(node.y),:);
node_SE = node(node.x == x_max & node.z == min(node.z) & node.y == max(node.y),:);

% base shaded region
r_width = 1500; % Hard code to ICSB for now
r_height = 900;
x_center = (max(node.x)-min(node.x))/2 + min(node.x);
z_center = (max(node.z)-min(node.z))/2 + min(node.z);

% Find analysis twist coordinates
SW_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_SW.id) '.mat']);
NW_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_NW.id) '.mat']);
NE_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_NE.id) '.mat']);
SE_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_SE.id) '.mat']);

[val, peak_twist_idx] = max(abs(NE_TH.nd_TH.disp_z_TH - NW_TH.nd_TH.disp_z_TH));
rotation = val/r_width;

x_twist(1) = x_plt_scale*rotation*r_height/2 + node_SW.x;
x_twist(2) = -x_plt_scale*rotation*r_height/2 + node_NW.x;
x_twist(3) = -x_plt_scale*rotation*r_height/2 + node_NE.x;
x_twist(4) = x_plt_scale*rotation*r_height/2 + node_SE.x;
x_twist(5) = x_twist(1);
z_twist(1) = z_plt_scale*SW_TH.nd_TH.disp_z_TH(peak_twist_idx) + node_SW.z;
z_twist(2) = z_plt_scale*NW_TH.nd_TH.disp_z_TH(peak_twist_idx) + node_NW.z;
z_twist(3) = z_plt_scale*NE_TH.nd_TH.disp_z_TH(peak_twist_idx) + node_NE.z;
z_twist(4) = z_plt_scale*SE_TH.nd_TH.disp_z_TH(peak_twist_idx) + node_SE.z;
z_twist(5) = z_twist(1);

% Find Recorded Twist Coordinates
if ~isempty(record_edp)
    [val, peak_twist_idx] = max(abs(record_edp.disp_TH_roof_east.z - record_edp.disp_TH_roof_west.z));
    rotation = val/r_width;
    
    x_twist_rec(1) = x_plt_scale*rotation*r_height/2 + node_SW.x;
    x_twist_rec(2) = -x_plt_scale*rotation*r_height/2 + node_NW.x;
    x_twist_rec(3) = -x_plt_scale*rotation*r_height/2 + node_NE.x;
    x_twist_rec(4) = x_plt_scale*rotation*r_height/2 + node_SE.x;
    x_twist_rec(5) = x_twist_rec(1);
    z_twist_rec(1) = z_plt_scale*record_edp.disp_TH_roof_west.z(peak_twist_idx) + node_SW.z;
    z_twist_rec(2) = z_plt_scale*record_edp.disp_TH_roof_west.z(peak_twist_idx) + node_NW.z;
    z_twist_rec(3) = z_plt_scale*record_edp.disp_TH_roof_east.z(peak_twist_idx) + node_NE.z;
    z_twist_rec(4) = z_plt_scale*record_edp.disp_TH_roof_east.z(peak_twist_idx) + node_SE.z;
    z_twist_rec(5) = z_twist_rec(1);
end

%% Plot
hold on
rectangle('Position',[x_min,min(node.z),r_width,r_height],'FaceColor',[0.85 0.85 0.85],'EdgeColor',[1 1 1])
set(gca,'XTick',[],'YTick',[])
plot_width = max([r_height,r_width])*1.25;
xlim([x_center-plot_width/2,x_center+plot_width/2])
ylim([z_center-plot_width/2,z_center+plot_width/2])

% Plot max twist from analaysis
plot(x_twist,z_twist,'LineWidth',2,'DisplayName','Analysis')

% Plot max twist from recording
if ~isempty(record_edp)
    plot(x_twist_rec,z_twist_rec,'--k','LineWidth',2,'DisplayName','Recorded')
end

% Format and save plot
axis off
h = legend;
pos = get(h,'Position');
posx = 0.45;
posy = 0.5;
set(h,'Position',[posx posy pos(3) pos(4)]);
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end