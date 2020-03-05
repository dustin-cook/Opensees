function [status] = driver_run_IDA(model_name, analysis_name, element_ids, node_ids, primary_nodes, story_ht, period)
% Script to Run and IDA of a building with a single ground motion %%


% Assumptions
% 1) Assumes gravity loads are defined and run in the model file
% 2) also assumes damping is defined in the model file

%% User Inputs
% Change user inputs to correct structure
try
    element = str2double(strsplit(element_ids,','));
    node.id = str2double(strsplit(node_ids,','));
    node.primary_nodes = str2double(strsplit(primary_nodes,','));
    story.story_ht = str2double(strsplit(story_ht,','));
    ida_results.period = str2double(strsplit(period,','));
catch
    status = 'IDA Failed: Issue with user inputs';
    return
end

% Ground Motion
analysis.gm_set = 'FEMA_far_field';

% Analysis options
analysis.summit = 0;
analysis.run_parallel = 0;
analysis.scale_increment = 0.25;
analysis.collapse_drift = 0.10;
analysis.clear_existing_data = 0;
analysis.run_ida = 1;
analysis.post_process_ida = 1;
analysis.general_ida = 1;
analysis.plot_ida = 1;

% Secondary options
analysis.opensees_SP = 1;
analysis.type = 1;
analysis.run_eigen = 0;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 0;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';

% Things that should be user inputs but are hardcoded
model.dimension = '2D';
analysis.run_z_motion = 0;
analysis.g_unit = 9810;

% Things to make fit with current framework
node.primary_story = ismember(node.id, node.primary_nodes);
hinge = [];
joint = [];

%% Import Packages
import ida.*

%% Define read and write directories
main_dir = ['outputs' '/' model_name '/' analysis_name];
tcl_dir = [main_dir '/' 'opensees_data'];

%% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);

%% Run Opensees Models
try
    fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
catch
    status = 'IDA Failed: IDA could not run';
    return
end

%% Collect IDA data and plot
try
    write_dir = [main_dir '/' 'IDA' '/' 'Summary Data'];
    
    % Create IDA table and plot IDA plot and collapse fragility
    if analysis.plot_ida
        fn_collect_and_plot_ida_simple(analysis, gm_set_table, node, story, main_dir, write_dir)
    end
    
    % Compresss summary folder
    zip([main_dir filesep 'IDA' filesep 'summary_data'],write_dir)
catch
    status = 'IDA Failed: IDA could not postprocess summary results';
    return
end

% If you get to here it means process completed successfully
status = 'IDA Completed Successfully';

end
