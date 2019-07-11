function [ ] = main_write_tcl( model_dimension, write_dir, node, element, story, joint, hinge, analysis, read_dir_analysis, ground_motion, hinge_grouping )
% Function to write TCL files for Opensees

%% Initial Setup
import opensees.write_tcl.*

%% Build Scripts
[ joint_ele_ids ] = fn_define_model( write_dir, node, element, joint, hinge, analysis, model_dimension, story, read_dir_analysis );
fn_define_recorders( write_dir, model_dimension, node, element, joint, hinge, analysis );
fn_define_loads( write_dir, analysis, node, model_dimension, story, element, ground_motion);
primary_nodes = node.id(node.primary_story == 1);
if analysis.run_eigen
    fn_eigen_analysis( write_dir, primary_nodes', length(story.id), analysis, model_dimension)
end
fn_setup_analysis( write_dir, write_dir, analysis, primary_nodes, story )
fn_define_analysis( write_dir, ground_motion, primary_nodes, story.story_ht, analysis, story )
end

