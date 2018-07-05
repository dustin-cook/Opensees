function [ ] = fn_plot_building( DCR, element, node, plot_name, plot_dir, dims, color, frame )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

if strcmp(frame,'ext')
    % Filter Elements to just be interior and exterior frame
    count = 0;
    for i = 1:length(element.id)
        if ~strcmp(element.type{i},'wall')
            if (node.z(node.id == element.node_1(i)) == 0) && (node.x(node.id == element.node_1(i)) >= 0) && (node.x(node.id == element.node_1(i)) <= 1500)
                count = count + 1;
                ele(count,:) = element(i,:);
                DCR_filter(count,:) = DCR(i,:);
            end
        end
    end

    node = node(node.z == 0 & node.x >=0 & node.x <= 1500,:);
elseif strcmp(frame,'int')
    % Filter Elements to just be interior and exterior frame
    count = 0;
    for i = 1:length(element.id)
        if ~strcmp(element.type{i},'wall')
            if (node.z(node.id == element.node_1(i)) == 300) && (node.x(node.id == element.node_1(i)) >= 0) && (node.x(node.id == element.node_1(i)) <= 1500)
                count = count + 1;
                ele(count,:) = element(i,:);
                DCR_filter(count,:) = DCR(i,:);
            end
        end
    end

    node = node(node.z == 300 & node.x >=0 & node.x <= 1500,:);
end

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

if strcmp(color,'raw')
    % Highlight elemets from 0 to 2 DCR
    s_break_1 = ele.node_1(DCR_filter < 2);
    t_break_1 = ele.node_2(DCR_filter < 2);
    highlight(H,s_break_1,t_break_1,'EdgeColor','g','LineWidth',1.5)

    % Highlight elemets from 2 to 4 DCR
    s_break_2 = ele.node_1(DCR_filter >= 2 & DCR_filter < 4);
    t_break_2 = ele.node_2(DCR_filter >= 2 & DCR_filter < 4);
    highlight(H,s_break_2,t_break_2,'EdgeColor','m','LineWidth',1.5)

    % Highlight elemets above 4 DCR
    s_break_3 = ele.node_1(DCR_filter >= 4);
    t_break_3 = ele.node_2(DCR_filter >= 4);
    highlight(H,s_break_3,t_break_3,'EdgeColor','r','LineWidth',1.5)

elseif strcmp(color,'linear')
    % Highlight elemets above 1 DCR acceptance
    s_break_3 = ele.node_1(DCR_filter >= 1);
    t_break_3 = ele.node_2(DCR_filter >= 1);
    highlight(H,s_break_3,t_break_3,'EdgeColor','r','LineWidth',2)
    
elseif strcmp(color,'nonlinear')
end

% Format and save plot
% if strcmp(dims,'3D')
%     xlabel('E-W Direcition (ft)')
%     ylabel('N-S Direcition (ft)')
%     zlabel('Height (ft)')
% else
    xlabel('Base (ft)')
    ylabel('Height (ft)')
% end
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end

