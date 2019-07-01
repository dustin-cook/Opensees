function [ ] = fn_plot_building_nl_3d( hinge, element, node, elev_title, plot_dir, direction, x_start, x_end, z_start, z_end )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp(direction,'x')
    type_filter = ~strcmp(element.type,'wall');
elseif strcmp(direction,'z')
    type_filter = strcmp(element.type,'wall');
end
        
% Filter Elements
count = 0;
for i = 1:length(element.id)
    if (node.z(node.id == element.node_1(i)) >= z_start) && (node.z(node.id == element.node_1(i)) <= z_end) &&(node.x(node.id == element.node_1(i)) >= x_start) && (node.x(node.id == element.node_1(i)) <= x_end)
        if (node.z(node.id == element.node_2(i)) >= z_start) && (node.z(node.id == element.node_2(i)) <= z_end) && (node.x(node.id == element.node_2(i)) >= x_start) && (node.x(node.id == element.node_2(i)) <= x_end)
            if type_filter(i)
                count = count + 1;
                ele(count,:) = element(i,:);
                ele_new.id(count) = ele.id(count);
                ele_new.old_node_1(count) = ele.node_1(count);
                ele_new.old_node_2(count) = ele.node_2(count);
                ele_new.node_1(count) = count*2-1;
                ele_new.node_2(count) = count*2;
                new_node.id(count*2-1) = count*2-1;
                new_node.id(count*2) = count*2;
                new_node.x(count*2-1) = node.x(node.id == ele.node_1(count));
                new_node.x(count*2) = node.x(node.id == ele.node_2(count));
                new_node.y(count*2-1) = node.y(node.id == ele.node_1(count));
                new_node.y(count*2) = node.y(node.id == ele.node_2(count));
                new_node.z(count*2-1) = node.z(node.id == ele.node_1(count));
                new_node.z(count*2) = node.z(node.id == ele.node_2(count));
            end
        end
    end
end

% Define Colors
cmap1 = [200 200 200;
        238 235 112;
        252 160 98;
        241 95 95]/255;
    
cmap2 = [0,60,48;
      1,102,94;
      53,151,143;
      128,205,193;
      199,234,229;
      246,232,195;
      223,194,125;
      191,129,45;
      140,81,10;
      84,48,5]/255;
    
%% Plot DCR on elevation of building
plot_elevation(ele_new, direction, new_node, hinge, cmap1, 4, plot_dir, ['ASCE 41 Acceptance - ' elev_title], 'accept')
plot_elevation(ele_new, direction, new_node, hinge, cmap2, 1, plot_dir, ['b Hinge Rotation - ' elev_title], 'b_ratio')



function [] = plot_elevation(ele_new, direction, new_node, hinge, cmap, cmap_max, plot_dir, plot_name, target)
% Import Packages
import plotting_tools.*

% Graph structure
s = ele_new.node_1;
t = ele_new.node_2;
G = graph(s,t);

% Plot Graph
if strcmp(direction,'z')
    x = new_node.z;
elseif strcmp(direction,'x')
    x = new_node.x;
end
y = new_node.y;
H = plot(G,'XData',x,'YData',y,'NodeLabel',{});
axis off

%% Highlight Performance
% Highlight elemets that pass IO
for e = 1:length(ele_new.id)
    hinges = hinge(hinge.element_id == ele_new.id(e) & strcmp(hinge.direction,'primary'),:);
    for h = 1:height(hinges)
        hinge_new_node_1 = ele_new.node_1(ele_new.old_node_1 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
        hinge_new_node_2 = ele_new.node_2(ele_new.old_node_2 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
        hinge_new_node = max([hinge_new_node_1, hinge_new_node_2]);
        cmap_length = length(cmap(:,1));
        cmap_idx = max([min([round(hinges.(target)(h)*cmap_length/cmap_max),cmap_length]),1]);
        highlight(H,hinge_new_node) 
        highlight(H,hinge_new_node,'NodeColor',cmap(cmap_idx,:)) % Highlight elemets that fail IO

    end
end

%% Format and save plot
title(plot_name)
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end
end

