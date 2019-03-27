function [ ] = fn_plot_building_nl_3d( hinge, element, node, plot_name, plot_dir, direction, x_start, x_end, z_start, z_end )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.*

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

%% Plot DCR on elevation of building
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

%% Highlight Performance
% Highlight elemets that pass IO
for e = 1:length(ele_new.id)
    hinges = hinge(hinge.element_id == ele_new.id(e) & strcmp(hinge.direction,'primary'),:);
    for h = 1:height(hinges)
        hinge_new_node_1 = ele_new.node_1(ele_new.old_node_1 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
        hinge_new_node_2 = ele_new.node_2(ele_new.old_node_2 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
        hinge_new_node = max([hinge_new_node_1, hinge_new_node_2]);
        
        if hinges.accept(h) == 2
            highlight(H,hinge_new_node) 
            highlight(H,hinge_new_node,'NodeColor','g') % Highlight elemets that fail IO
        elseif hinges.accept(h) == 3
            highlight(H,hinge_new_node)
            highlight(H,hinge_new_node,'NodeColor','y') % Highlight elemets that fail LS
        elseif hinges.accept(h) == 4
            highlight(H,hinge_new_node)
            highlight(H,hinge_new_node,'NodeColor','r') % Highlight elemets that fail CP
        end
    end
end


% element_list = hinge.element_id(hinge.accept == 1);
% for i = 1:length(element_list)
%     if (sum(ele.id == element_list(i))>0)
%         s_break_3 = ele_new.node_1(ele.id == element_list(i));
%         t_break_3 = ele_new.node_2(ele.id == element_list(i));
%         highlight(H,s_break_3,t_break_3,'EdgeColor','g','LineWidth',2)
%     end
% end
% 
% % Highlight elemets that pass LS
% element_list = hinge.element_id(hinge.accept == 2);
% for i = 1:length(element_list)
%     if (sum(ele.id == element_list(i))>0)
%         s_break_3 = ele_new.node_1(ele.id == element_list(i));
%         t_break_3 = ele_new.node_2(ele.id == element_list(i));
%         highlight(H,s_break_3,t_break_3,'EdgeColor','b','LineWidth',2)
%     end
% end
% 
% % Highlight elemets that pass CP
% element_list = hinge.element_id(hinge.accept == 3);
% for i = 1:length(element_list)
%     if (sum(ele.id == element_list(i))>0)
%         s_break_3 = ele_new.node_1(ele.id == element_list(i));
%         t_break_3 = ele_new.node_2(ele.id == element_list(i));
%         highlight(H,s_break_3,t_break_3,'EdgeColor','y','LineWidth',2)
%     end
% end
% 
% % Highlight elemets that fail all
% element_list = hinge.element_id(hinge.accept == 4);
% for i = 1:length(element_list)
%     if (sum(ele.id == element_list(i))>0)
%         s_break_3 = ele_new.node_1(ele.id == element_list(i));
%         t_break_3 = ele_new.node_2(ele.id == element_list(i));
%         highlight(H,s_break_3,t_break_3,'EdgeColor','r','LineWidth',2)
%     end
% end

%% Format and save plot
xlabel('Base (ft)')
ylabel('Height (ft)')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end

