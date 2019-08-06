function [] = fn_run_gm_stripe(analysis, model, story, element, node, hinge, scale_factor, gm_set_table, gms, ida_results, tcl_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_main_IDA
import ida.fn_postprocess_ida

% Run Opensees
if analysis.run_ida
    fprintf('Running Scale Factor %4.2f for Ground Motion ID: %i-%i \n', scale_factor, gm_set_table.set_id(gms), gm_set_table.pair(gms))
    [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor, ida_results.period, tcl_dir);
    if exit_status == 1
        printf('Failed GMs set ID: %i and Pair: %i for Scale Factor %4.2f \n', gm_set_table.set_id(gms), gm_set_table.pair(gms), scale_factor)
    end
else
    exit_status = 0;
end

if analysis.post_process_ida && exit_status ~= 1
    fprintf('Postprocessing Opensees Ouputs\n')
    fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gms, scale_factor)
end
end

