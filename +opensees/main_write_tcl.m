function [ node, ground_motion ] = main_write_tcl( model, output_dir, node, element, story, joint, hinge, analysis )
% Function to write TCL files for Opensees

%% Initial Setup
import opensees.write_tcl.*

%% Build Scripts
[ node ] = fn_define_model( output_dir, node, element, joint, hinge, analysis, model.dimension, story );
fn_define_recorders( output_dir, model.dimension, node.id', element.id, hinge, analysis )
[ground_motion] = fn_define_loads( output_dir, analysis, node, model.dimension, story, element.id');
first_story_node = node.id(node.primary_story == 1);
if analysis.run_eigen
    fn_eigen_analysis( output_dir, analysis.time_step, first_story_node', length(story.id), analysis)
end
if analysis.type == 1 % Dynamic
    fn_setup_dynamic_analysis( output_dir, analysis )
elseif analysis.type == 2 % Pushover
    fn_setup_pushover_analysis( output_dir, analysis, node )
end
fn_define_analysis( output_dir, ground_motion, first_story_node, story.story_ht, analysis )
end

