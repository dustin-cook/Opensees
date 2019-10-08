function [ ] = fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_run_gm_ida

if analysis.run_parallel
    parpool; % Set up Parallel Workers
end

% Loop through each ground motion
tim_start = tic;
if analysis.run_parallel
    parfor gms = 1:height(gm_set_table)
        % Loop through each scale of the GM
        fn_run_gm_ida(analysis, model, story, element, node, hinge, joint, gm_set_table, gms, ida_results, tcl_dir)
    end
else
    for gms = 1:height(gm_set_table)
        % Loop through each scale of the GM
        fn_run_gm_ida(analysis, model, story, element, node, hinge, joint, gm_set_table, gms, ida_results, tcl_dir)
    end
end
tim_elapsed = toc(tim_start);
fprintf('IDA finished with a run time of %4.2f seconds \n', tim_elapsed)

delete(gcp('nocreate')) % End Any Parallel Process
end

