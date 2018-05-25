function [ ] = fn_plot_building( DCR, element, node, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Graph structure
s = element.node_1;
t = element.node_2;
G = graph(s,t);

% Reduce nodes to only ones that matter for elements 
for i = 1:max([element.node_1;element.node_2])
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

% Highlight elemets from 0.75 to 1 DCR
s_break_1 = element.node_1(DCR >= 0.75 & DCR < 1);
t_break_1 = element.node_2(DCR >= 0.75 & DCR < 1);
highlight(H,s_break_1,t_break_1,'EdgeColor','g','LineWidth',1.5)

% Highlight elemets from 1 to 2 DCR
s_break_2 = element.node_1(DCR >= 1 & DCR < 2);
t_break_2 = element.node_2(DCR >= 1 & DCR < 2);
highlight(H,s_break_2,t_break_2,'EdgeColor','m','LineWidth',2)

% Highlight elemets above 2 DCR
s_break_3 = element.node_1(DCR >= 2);
t_break_3 = element.node_2(DCR >= 2);
highlight(H,s_break_3,t_break_3,'EdgeColor','r','LineWidth',2.5)

% Format and save plot
xlabel('Base (ft)')
ylabel('Height (ft)')
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end

