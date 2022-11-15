function [ ] = fn_curt_plot_2D( hinge, element, node, story, plot_name, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import plotting_tools.*

% Filter Elements
count = 0;
for i = 1:length(element.id)
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
    % determine displacements
%     if strcmp(ele.type(count),'column') || strcmp(ele.type(count),'wall')
%         story_disp = max([story.max_disp_x(story.id == (ele.story(count)-1)),0]);
%     else
%         story_disp = max([story.max_disp_x(story.id == ele.story(count)),0]);
%     end
    story_disp = max([story.max_disp_x(story.id == ele.story(count)),0]);
    story_below_disp = max([story.max_disp_x(story.id == (ele.story(count)-1)),0]);
    rel_story_disp = story_disp - story_below_disp;
    story_y_top = max([story.y_start(story.id == (ele.story(count))) + story.story_ht(story.id == (ele.story(count))),0]);
    story_y_bot = max([story.y_start(story.id == (ele.story(count)-1)) + story.story_ht(story.id == (ele.story(count)-1)),0]);
    story_y_rel = story_y_top - story_y_bot;
    node_y_rel_1 = new_node.y(count*2-1) - story_y_bot;
    node_y_rel_2 = new_node.y(count*2) - story_y_bot;
    rel_node_disp_1 = rel_story_disp*(node_y_rel_1/story_y_rel);
    rel_node_disp_2 = rel_story_disp*(node_y_rel_2/story_y_rel);
    new_node.x_disp(count*2-1) = story_below_disp + rel_node_disp_1;
    new_node.x_disp(count*2) = story_below_disp + rel_node_disp_2;
end

% Graph structure
s = ele_new.node_1;
t = ele_new.node_2;
G = graph(s,t);

% Plot Graph
H = plot(G,'XData',new_node.x+20*new_node.x_disp,'YData',new_node.y,'NodeLabel',{});
axis off

%% Highlight Performance
cmap = jet(100);

% Highlight elemets that pass Acceptance Criteria
for e = 1:length(ele_new.id)
    hinges = hinge(hinge.element_id == ele_new.id(e) & strcmp(hinge.direction,'primary'),:);
    for h = 1:height(hinges)
        if any(strcmp('node',hinges.Properties.VariableNames))
            hinge_node_1 = ele_new.node_1(ele_new.old_node_1 == hinges.node(h));
            hinge_node_2 = ele_new.node_2(ele_new.old_node_2 == hinges.node(h));
            hinge_node = max([hinge_node_1, hinge_node_2]);
        else
            hinge_node_1 = ele_new.node_1(ele_new.old_node_1 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
            hinge_node_2 = ele_new.node_2(ele_new.old_node_2 == hinges.node_1(h) | ele_new.old_node_1 == hinges.node_2(h));
            hinge_node = max([hinge_node_1, hinge_node_2]);
        end
%         if hinges.a_ratio(h) >= 1
%             highlight(H,hinge_node,'MarkerSize',min([max([hinges.b_ratio(h)*5,2]),10])) % Make hinge bigger
%             highlight(H,hinge_node,'NodeColor','r') % Highlight hinges
%         end
        
        % Just highlight based on a-ratio
        c_idx = min(ceil(hinges.plastic_deform(h) * 10000)+1,100);
        highlight(H,hinge_node,'NodeColor', cmap(c_idx,:)) % Highlight hinges
    end
end

%% Format and save plot
fn_format_and_save_plot( plot_dir, plot_name, 4 )

end

