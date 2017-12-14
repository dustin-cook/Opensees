function [ node, element ] = fn_model_table( num_stories, num_bays, floor_ht, bay_width, foundation_fix, story_mass, story_weight, story_force, A_col, E_col, I_col, A_bm, E_bm, I_bm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Create Governing Parameters
num_nodes = (num_bays + 1) * (num_stories + 1);
nodes = 1:num_nodes;
nodal_coordinates_x = (0:num_bays)*bay_width;
nodal_coordinates_y = floor_ht;

% Create Coordinate Vectors
nodal_y = [];
nodal_x = [];
for i = 1:length(nodal_coordinates_y)
    nodal_x = [nodal_x,nodal_coordinates_x];
    nodal_y = [nodal_y,nodal_coordinates_y(i)*ones(1,length(nodal_coordinates_x))];
end

% Define Nodal Fixity, Mass, Weight, and force
node_id = 0;
for story = 0:num_stories
    for bay = 1:num_bays+1
        node_id = node_id+1;
        if story == 0 % Ground Nodes
            nodal_fix{node_id} = foundation_fix;
            nodal_mass(node_id) = 0;
            nodal_wt(node_id) = 0;
            nodal_force(node_id) = 0;
        elseif bay == 1 || bay == num_bays+1 % Frame end nodes
            nodal_fix{node_id} = [0,0,0];
            nodal_mass(node_id) = story_mass(story)/(num_bays*2);
            nodal_wt(node_id) = story_weight(story)/(num_bays*2);
            nodal_force(node_id) = story_force(story)/(num_bays*2);
        else
            nodal_fix{node_id} = [0,0,0];
            nodal_mass(node_id) = story_mass(story)/(num_bays);
            nodal_wt(node_id) = story_weight(story)/(num_bays);
            nodal_force(node_id) = story_force(story)/(num_bays);
        end
    end
end

% Define Elements 
num_beam = num_bays*num_stories;
num_col = (num_bays+1)*num_stories;
num_ele = num_beam + num_col;
elements = 1:num_ele;

col_node_start = 1:(num_nodes-(num_bays+1));
col_node_end = (1+(num_bays+1)):num_nodes;
logic_points = ones(num_bays+1,num_stories+1);
logic_points(:,1) = 0;
logic_points(end,:) = 0;
beam_node_start = find(logic_points)';

logic_points = ones(num_bays+1,num_stories+1);
logic_points(:,1) = 0;
logic_points(1,:) = 0;
beam_node_end = find(logic_points)';

node_start = [col_node_start,beam_node_start];
node_end = [col_node_end,beam_node_end];

% Define Element Properties
area_col = [];
youngs_mod_col = [];
moment_of_inertia_col = [];
area_bm = [];
youngs_mod_bm = [];
moment_of_inertia_bm = [];
for n = 1:num_stories
    area_col = [area_col,ones(1,num_bays+1)*A_col(n)];
    youngs_mod_col = [youngs_mod_col,ones(1,num_bays+1)*E_col(n)];
    moment_of_inertia_col = [moment_of_inertia_col,ones(1,num_bays+1)*I_col(n)];
    area_bm = [area_bm,ones(1,num_bays)*A_bm(n)];
    youngs_mod_bm = [youngs_mod_bm,ones(1,num_bays)*E_bm(n)];
    moment_of_inertia_bm = [moment_of_inertia_bm,ones(1,num_bays)*I_bm(n)];
end
area = [area_col,area_bm];
youngs_mod = [youngs_mod_col,youngs_mod_bm];
moment_of_inertia = [moment_of_inertia_col,moment_of_inertia_bm];

% Create Model Data Table
node = table(nodes',nodal_x',nodal_y',nodal_fix',nodal_mass',nodal_wt',nodal_force','VariableNames',{'id','x_coor','y_coor','fix','mass','weight','force'});
element = table(elements',node_start',node_end',area',youngs_mod',moment_of_inertia','VariableNames',{'id','node_start','node_end','A','E','I'});

end

