function [ ] = fn_master_IDA(analysis, model, story, element, node, hinge, IDA_scale_factors,gm_set_table, ida_results, tcl_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_main_IDA
import ida.fn_postprocess_ida

parpool; % Set up Parallel Workers
for i = 1:length(IDA_scale_factors)
    error_count = 0;
    scale_factor = IDA_scale_factors(i);
    analysis.ground_motion_scale_factor = scale_factor;
    run_ida = analysis.run_ida;
    parfor gms = 1:height(gm_set_table)
        % Run Opensees
        if run_ida
            % Suppress MATLAB warnings
            warning('off','all')
            fprintf('Running Scale Factor %4.2f for Ground Motion ID: %i-%i \n\n', scale_factor, gm_set_table.set_id(gms), gm_set_table.pair(gms))
            [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, ida_results.period, tcl_dir);
            if exit_status == 1
                error_count = error_count + 1;
            end
        else
            exit_status = 0;
        end

        if analysis.post_process_ida && exit_status ~= 1
            fprintf('Postprocessing Opensees Ouputs\n')
            fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor)
        end
        fprintf('\n')
    end
    fprintf('%i Failed GMs for Scale Factor %4.2f \n\n', error_count, scale_factor)
end
delete(gcp('nocreate')) % End Parallel Process
end

