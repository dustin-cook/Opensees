function [ node, element, story_ht, num_node_story ] = fn_model_table( model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
story_force = zeros(1,model.num_stories); % Just placeholder for now, need to fix later

%% Load data tables
story_table = readtable(['inputs' filesep 'story.csv'],'ReadVariableNames',true);
element_group_table = readtable(['inputs' filesep 'element_group.csv'],'ReadVariableNames',true);
element_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
node_group_table = readtable(['inputs' filesep 'node_group.csv'],'ReadVariableNames',true);

% Object Methods
story = story_table(story_table.model == model.id,:);
floor_ht = zeros(model.num_stories,1);
for i = 1:model.num_stories
    story_ht(i) = story.ht(story.story_index == i,:);
    floor_ht(i+1) = story_ht(i) + floor_ht(i);
    story_wt(i) = story.wt(story.story_index == i,:);
end
story_mass = story_wt/386;
building_ht = sum(story_ht);
building_mass = sum(story_mass);
if strcmp(model.foundation,'fix')
    foundation_fix = [1 1 1];
elseif strcmp(model.foundation,'pin')
    foundation_fix = [1 1 0];
elseif strcmp(model.foundation,'roller')
    foundation_fix = [0 1 0];
else % Assume fixed foundation support
    foundation_fix = [1 1 1];
end

%% Loop through Stories to create nodes and elements
num_node_story = model.num_bays*5 + 3;  % Number of nodes in each story
num_bms_story = model.num_bays*3;
num_col_story_typ = model.num_bays*3 + 3;
num_col_story_ground = model.num_bays*2 + 2;
num_nodes = num_node_story * model.num_stories;
num_bms = num_bms_story * model.num_stories;
num_col = num_col_story_typ * (model.num_stories-1) + num_col_story_ground;
nodal_x = zeros(1,num_nodes);
nodal_y = zeros(1,num_nodes);
nodal_fix = cell(1,num_nodes);
nodal_mass = zeros(1,num_nodes);
nodal_wt = zeros(1,num_nodes);
nodal_force = zeros(1,num_nodes);
bm_a_array = [];
bm_e_array = [];
bm_i_array = [];
bm_node_start = [];
bm_node_end = [];
col_a_array = [];
col_e_array = [];
col_i_array = [];
col_node_start = [];
col_node_end = [];

for s = 1:model.num_stories
    % Define Critical Nodes
    first_node = (s-1) * num_node_story + 1;
    layer_2 = first_node + model.num_bays + 1;
    layer_3 = layer_2 + model.num_bays + 1;
    last_node = s * num_node_story;
    layer_3_center_array = layer_3:3:last_node;
    
    % Filter to get Element Grouping for this Story
    element_group = element_group_table(element_group_table.id == story.element_group_id(s),:);
    node_group = node_group_table(node_group_table.id == story.node_group_id(s),:);
    
    % Assemble x coordinate of nodes based on joint size
    x_coor_story_array = [];
    bm_d_story_array_upper = [];
    bm_d_story_array_lower = [];

    for i = 1:model.num_bays
        col_d_1 = element_table.d(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 90));
        col_d_2 = element_table.d(element_table.id == element_group.element(element_group.position == (i+1) & element_group.orientation == 90));
        bm_d = element_table.d(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 0));
        if i == 1
            bm_d_story_array_upper = [bm_d_story_array_upper, bm_d, bm_d, bm_d, bm_d];
        else
            bm_d_story_array_upper = [bm_d_story_array_upper, bm_d, bm_d, bm_d];
        end
        bm_d_story_array_lower = [bm_d_story_array_lower, bm_d];
        x_coor_story_array = [x_coor_story_array, col_d_1/2, model.bay_width - (col_d_1/2+col_d_2/2), col_d_2/2];
        
        % Beam Element Assignment
        bm_a = element_table.a(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 0));
        bm_e = element_table.e(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 0));
        bm_i = element_table.i(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 0));
        bm_a_array = [bm_a_array, 9e6, bm_a, 9e6];
        bm_e_array = [bm_e_array, bm_e, bm_e, bm_e]; 
        bm_i_array = [bm_i_array, 9e9, bm_i, 9e9];
    end
    bm_node_start = [bm_node_start, (layer_3):(last_node-1)];
    bm_node_end = [bm_node_end, (layer_3+1):(last_node)];
    
    for i = 1:model.num_bays+1
        % Column Element Assignment
        col_a(i) = element_table.a(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 90));
        col_e(i) = element_table.e(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 90));
        col_i(i) = element_table.i(element_table.id == element_group.element(element_group.position == i & element_group.orientation == 90));
    end
    if s == 1
        col_a_array = [col_a_array, col_a, col_a];
        col_e_array = [col_e_array, col_e, 100*col_e]; 
        col_i_array = [col_i_array, col_i, col_i];
        col_node_start = [col_node_start, (first_node):(layer_2-1), (layer_2):(layer_3-1)];
        col_node_end = [col_node_end, (layer_2):(layer_3-1), layer_3_center_array];
    else
        col_a_array = [col_a_array, col_a, col_a, col_a];
        col_e_array = [col_e_array, 100*col_e, col_e, 100*col_e]; 
        col_i_array = [col_i_array, col_i, col_i, col_i];
        col_node_start = [col_node_start, col_node_end((end-model.num_bays):end), (first_node):(layer_2-1), (layer_2):(layer_3-1)];
        col_node_end = [col_node_end, (first_node):(layer_2-1), (layer_2):(layer_3-1), layer_3_center_array];
    end
    
    % Assign Nodal Properties
    for n = first_node:last_node
        % X and Y Coordinate
        if n >= layer_3
            if model.num_bays == 0
                nodal_y(n) = sum(story_ht(1:s)) - 1/2;
            else
                nodal_y(n) = sum(story_ht(1:s)) - bm_d_story_array_upper(n-(layer_3-1))/2;
            end
            if n == layer_3
                nodal_x(n) = 0;
            else
                nodal_x(n) = nodal_x(n-1)+x_coor_story_array(n-layer_3);
            end
        elseif n >= layer_2
            nodal_x(n) = (n-layer_2) * model.bay_width;
            if n == layer_2 && model.num_bays == 0
                nodal_y(n) = sum(story_ht(1:s)) - 1;
            elseif n == layer_2
                nodal_y(n) = sum(story_ht(1:s)) - bm_d_story_array_lower(1);
            elseif n == (layer_3-1)
                nodal_y(n) = sum(story_ht(1:s)) - bm_d_story_array_lower(end);
            else
                nodal_y(n) = sum(story_ht(1:s)) - max(bm_d_story_array_lower(n-(layer_2)),bm_d_story_array_lower(n-(layer_2-1)));
                
            end
        else
            nodal_y(n) = sum(story_ht(1:(s-1)));
            nodal_x(n) = (n-first_node) * model.bay_width; 
        end
        
        % Foundation Fixity
        if s == 1 && n < layer_2
            nodal_fix{n} = foundation_fix;
        else
            nodal_fix{n} = [0, 0, 0];
        end
        
        % Nodal Mass and Weight and Applied Nodal Force
        if n >= layer_3 && sum(layer_3_center_array == n) > 0
            if model.num_bays == 0
                nodal_mass(n) = story_mass(s);
                nodal_wt(n) = story_wt(s);
                nodal_force(n) = story_force(s);
            else
                trib_ratio = node_group.trib_ratio(node_group.node == (((n-layer_3)/3)+1));
                nodal_mass(n) = (story_mass(s)/(model.num_bays))*trib_ratio;
                nodal_wt(n) = (story_wt(s)/(model.num_bays))*trib_ratio;
                nodal_force(n) = (story_force(s)/(model.num_bays))*trib_ratio;
            end
        else
            nodal_mass(n) = 0;
            nodal_wt(n) = 0;
            nodal_force(n) = 0;
        end
    end
end

% Prep assigments for table assembly
nodes = 1:num_nodes;
elements = 1:(num_col+num_bms);
area = [col_a_array,bm_a_array];
youngs_mod = [col_e_array,bm_e_array];
moment_of_inertia = [col_i_array,bm_i_array];
node_start = [col_node_start,bm_node_start];
node_end = [col_node_end,bm_node_end];

% Create Model Data Table
node = table(nodes',nodal_x',nodal_y',nodal_fix',nodal_mass',nodal_wt',nodal_force','VariableNames',{'id','x_coor','y_coor','fix','mass','weight','force'});
element = table(elements',node_start',node_end',area',youngs_mod',moment_of_inertia','VariableNames',{'id','node_start','node_end','A','E','I'});

end

