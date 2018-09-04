function [ node, ground_motion ] = main_write_tcl( model, output_dir, node, element, story, joint, hinge, analysis )
% Function to write TCL files for Opensees

%% Initial Setup
import opensees.write_tcl.*

%% Build Scripts
if strcmp(model.dimension,'2D')
    [ node ] = fn_build_model_2D( output_dir, node, element, joint, hinge, analysis, model.dimension );
elseif strcmp(model.dimension,'3D')
    [ node ] = fn_build_model_3D( output_dir, node, element, story, joint, hinge, analysis );
else
    error('Number of Dimensions Not Recognized')
end
fn_define_recorders( output_dir, model.dimension, node.id', element.id, hinge, analysis )
[ground_motion] = fn_define_loads( output_dir, analysis, node, model.dimension, length(story.id), element.id');
first_story_node = node.id(node.primary_story == 1);
if analysis.run_eigen
    fn_eigen_analysis( output_dir, analysis.time_step, first_story_node', length(story.id), analysis)
end
fn_setup_dynamic_analysis( output_dir, analysis )
fn_define_analysis( output_dir, ground_motion, first_story_node, story.story_ht, analysis )
end

