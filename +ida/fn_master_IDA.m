function [ ] = fn_master_IDA(analysis, model, story, element, node, hinge, IDA_scale_factors,gm_set_table, ida_results, tcl_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_run_gm_stripe

if analysis.run_parallel
    parpool; % Set up Parallel Workers
end

% Loop through each stripe
for i = 1:length(IDA_scale_factors)
    scale_factor = IDA_scale_factors(i);
    analysis.ground_motion_scale_factor = scale_factor;
    
    % Loop through each ground motion in the stripe
    if analysis.run_parallel
        parfor gms = 1:height(gm_set_table)
            fn_run_gm_stripe(analysis, model, story, element, node, hinge, scale_factor, gm_set_table, gms, ida_results, tcl_dir)
        end
    else
        for gms = 1:height(gm_set_table)
            fn_run_gm_stripe(analysis, model, story, element, node, hinge, scale_factor, gm_set_table, gms, ida_results, tcl_dir)
        end
    end

end

delete(gcp('nocreate')) % End Any Parallel Process
end

