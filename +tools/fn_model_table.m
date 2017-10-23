function [ node, element ] = fn_model_table( num_stories, num_bays, story_ht, bay_width, foundation_fix, story_mass, story_weight, story_force, A, E, I )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Create Governing Parameters
num_nodes = (num_bays + 1) * (num_stories + 1);
nodes = 1:num_nodes;
nodal_coordinates_x = (0:num_bays)*bay_width;
nodal_coordinates_y = (0:num_stories)*story_ht;

% Create Coordinate Vectors
nodal_y = [];
nodal_x = [];
for i = 1:length(nodal_coordinates_y)
    nodal_x = [nodal_x,nodal_coordinates_x];
    nodal_y = [nodal_y,nodal_coordinates_y(i)*ones(1,length(nodal_coordinates_x))];
end

% Define Nodal Fixity, Mass, Weight, and force
for i = 1:num_nodes
    if i <= (num_bays + 1)
        nodal_fix{i} = foundation_fix;
        nodal_mass(i) = 0;
    else
        nodal_fix{i} = [0,0,0];
        nodal_mass(i) = story_mass/(num_bays + 1);
        nodal_wt(i) = story_weight/(num_bays + 1);
        nodal_force(i) = story_force/(num_bays + 1);
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
area = ones(1,num_ele)*A;
youngs_mod = ones(1,num_ele)*E;
moment_of_inertia = ones(1,num_ele)*I;

% Create Model Data Table
node = table(nodes',nodal_x',nodal_y',nodal_fix',nodal_mass',nodal_wt',nodal_force','VariableNames',{'id','x_coor','y_coor','fix','mass','weight','force'});
element = table(elements',node_start',node_end',area',youngs_mod',moment_of_inertia','VariableNames',{'id','node_start','node_end','A','E','I'});

end

