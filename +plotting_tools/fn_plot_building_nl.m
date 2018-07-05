function [ ] = fn_plot_building_nl( hinge, ele, node, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

% Filter Elements to just be interior and exterior frame


%% Plot DCR on elevation of building
% Graph structure
s = ele.node_1;
t = ele.node_2;
G = graph(s,t);

% Reduce nodes to only ones that matter for elements 
% if strcmp(dims,'3D')
%     for i = 1:max([ele.node_1;ele.node_2])
%         if sum(node.id == i) == 0
%             x(i) = 0;
%             y(i) = 0;
%             z(i) = 0;
%         else
%             x(i) = node.x(node.id == i);
%             y(i) = node.y(node.id == i);
%             z(i) = node.z(node.id == i);
%         end
%     end
% 
%     % Plot Graph
%     H = plot(G,'XData',x,'YData',z,'ZData',y);
% else
    for i = 1:max([ele.node_1;ele.node_2])
        if sum(node.id == i) == 0
            x(i) = 0;
            y(i) = 0;
        else
            x(i) = node.x(node.id == i);
            y(i) = node.y(node.id == i);
        end
    end

    % Plot Graph
    H = plot(G,'XData',x,'YData',y);
% end

% Highlight elemets that pass IO
element_list = hinge.element_id(hinge.accept == 1);
for i = 1:length(element_list)
    s_break_3 = ele.node_1(ele.id == element_list(i));
    t_break_3 = ele.node_2(ele.id == element_list(i));
    highlight(H,s_break_3,t_break_3,'EdgeColor','b','LineWidth',2)
end

% Highlight elemets that pass LS
element_list = hinge.element_id(hinge.accept == 2);
for i = 1:length(element_list)
    s_break_3 = ele.node_1(ele.id == element_list(i));
    t_break_3 = ele.node_2(ele.id == element_list(i));
    highlight(H,s_break_3,t_break_3,'EdgeColor','g','LineWidth',2)
end

% Highlight elemets that pass CP
element_list = hinge.element_id(hinge.accept == 3);
for i = 1:length(element_list)
    s_break_3 = ele.node_1(ele.id == element_list(i));
    t_break_3 = ele.node_2(ele.id == element_list(i));
    highlight(H,s_break_3,t_break_3,'EdgeColor','y','LineWidth',2.5)
end


% Highlight elemets that fail all
element_list = hinge.element_id(hinge.accept == 4);
for i = 1:length(element_list)
    s_break_3 = ele.node_1(ele.id == element_list(i));
    t_break_3 = ele.node_2(ele.id == element_list(i));
    highlight(H,s_break_3,t_break_3,'EdgeColor','r','LineWidth',2)
end


% Format and save plot
xlabel('Base (ft)')
ylabel('Height (ft)')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end

