function [ ] = fn_build_sdof( model, output_dir )
% Funtion to build sdof element and node tables

%% INITIAL SETUP
% Import Packages
import build_model.fn_define_hinge

%% Build Story Table
story.id = 1;
story.story_ht = model.height;

%% Build Node Table
node.id = [1;2];
node.x = [0;0];
node.y = [0; model.height];
node.dead_load = [0; model.weight];
node.live_load = [0;0];
node.mass = node.dead_load/386;
node.fix = {'[111111]';'[000000]'};
node.story = [0;1];
node.primary_story = [0;1];

%% Build Element Table
element.id = 1;
element.ele_id = 0;
element.type = 'column';
element.node_1 = 1;
element.node_2 = 2;
element.length = model.height;
element.k = 2*pi*(model.weight/386)/(model.period^2);
element.e = 1000000;
element.a = 1000000;
element.i = element.k*(element.length^3)/(3*element.e);

%% Assign Joints
joint.id = [];

%% Create Nonlinear Rotational Springs at ends of all beams and columns
hinge.id = [];
hinge.element_id = [];
hinge.node_1 = [];
hinge.node_2 = [];
if model.nonlinear ~= 0
    [ node, element, hinge ] = fn_define_hinge( hinge, element, node, 0, 0 );
    element.Mn_aci_pos = model.moment_capacity;
    element.Mn_aci_neg = model.moment_capacity;
end

%% Reformat outputs to table and write CSV's
story_table = struct2table(story);
writetable(story_table,[output_dir filesep 'story.csv'])
node_table = struct2table(node);
writetable(node_table,[output_dir filesep 'node.csv'])
ele_table = struct2table(element);
writetable(ele_table,[output_dir filesep 'element.csv'])
joint_table = struct2table(joint);
writetable(joint_table,[output_dir filesep 'joint.csv'])
hinge_table = struct2table(hinge);
writetable(hinge_table,[output_dir filesep 'hinge.csv'])

end

