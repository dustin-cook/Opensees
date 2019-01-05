function [ node, ground_motion ] = main_write_tcl( model_dimension, write_dir, node, element, story, joint, hinge, analysis, read_dir_analysis )
% Function to write TCL files for Opensees

%% Initial Setup
import opensees.write_tcl.*

%% Build Scripts
[ node ] = fn_define_model( write_dir, node, element, joint, hinge, analysis, model_dimension, story, read_dir_analysis );
fn_define_recorders( write_dir, model_dimension, node.id', element, hinge, analysis )
[ground_motion] = fn_define_loads( write_dir, analysis, node, model_dimension, story, element.id');
first_story_node = node.id(node.primary_story == 1);
if analysis.run_eigen
    fn_eigen_analysis( write_dir, first_story_node', length(story.id), analysis)
end
if analysis.type == 1 || analysis.type == 2 % Dynamic or Pushover
    fn_setup_analysis( write_dir, analysis, first_story_node, story )
elseif analysis.type == 3 % Static Cyclic
    fn_setup_static_cyclic_analysis( write_dir, analysis )
end
fn_define_analysis( write_dir, ground_motion, first_story_node, story.story_ht, analysis, story )
end

