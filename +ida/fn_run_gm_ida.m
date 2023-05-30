function [] = fn_run_gm_ida(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gm_idx, ida_results, tcl_dir, main_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_main_IDA
import ida.fn_postprocess_ida
import ida.fn_postprocess_ida_gen

% Define scale factor array
scale_factor = 0;
scale_increment = analysis.scale_increment;
prev_collapse = 0;

% Run for each scale factor
while scale_factor < analysis.scale_increment*25
    % Increment Scale Factor
    scale_factor = scale_factor + scale_increment;
    analysis.ground_motion_scale_factor = scale_factor;
    
    % Define Ground Motion Data
    ground_motion.x = gms2run(gm_idx,:);
    ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
    ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
    if analysis.run_z_motion
        ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
        ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
        ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};
    end
    
    % Load spectral info and save Sa
    spectra_table = readtable([ground_motion.x.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
    ground_motion.x.sa = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1))*scale_factor;
    if analysis.run_z_motion
        spectra_table = readtable([ground_motion.z.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
        ground_motion.z.sa = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(2))*scale_factor;
    end

    % Create Output Directories
    ida_opensees_dir = [main_dir '/' 'IDA' '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair) '/' 'Scale_' num2str(scale_factor) ];
    if ~exist(ida_opensees_dir,'dir')
        mkdir(ida_opensees_dir)
    end

    ida_summary_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair) '/' 'Scale_' num2str(scale_factor)];
    if ~exist(ida_summary_dir,'dir')
        mkdir(ida_summary_dir)
    end

    % Run Opensees
    if analysis.run_ida
        fprintf('Running GM ID: %i Pair: %i for Scale Factor %4.2f \n', gms2run.set_id(gm_idx), gms2run.pair(gm_idx), scale_factor)
        tim_start = tic;
        [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, joint, ground_motion, tcl_dir, ida_opensees_dir, ida_summary_dir);
        tim_elapsed = toc(tim_start);
        if exit_status == 1
            fprintf('Failed GM ID: %i Pair: %i for Scale Factor %4.2f \n', gms2run.set_id(gm_idx), gms2run.pair(gm_idx), scale_factor)
        else
            fprintf('Successful run time of %4.2f seconds for GM ID: %i Pair: %i for Scale Factor %4.2f \n', tim_elapsed, gms2run.set_id(gm_idx), gms2run.pair(gm_idx), scale_factor)
        end
    else
        exit_status = 0;
    end

    % Post process single IDA run
    if analysis.post_process_ida && exit_status ~= 1
        fprintf('Postprocessing Opensees Ouputs\n')
        if analysis.general_ida
            fn_postprocess_ida_gen( analysis, model, node, element, ground_motion, ida_opensees_dir, ida_summary_dir )
        else 
            fn_postprocess_ida(analysis, model, story, element, node, hinge, joint, ground_motion, ida_opensees_dir, ida_summary_dir, ele_props_table)
        end
    end
    
    % Check if collapse was reached
    load([ida_summary_dir filesep 'summary_results.mat'])
    if summary.collapse > 0
        scale_increment = -abs(scale_increment/2);
        prev_collapse = 1;
    elseif  prev_collapse
        scale_increment = abs(scale_increment/2);
    end
    
    % Check if scale increment is small enough to be negligable
    if abs(scale_increment) < 0.05
        break
    end
end

% Write flag idicating that the ground motion has finished
fileID = fopen([main_dir '/' 'IDA' '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair) '/' 'gm_complete.txt'],'w');
fprintf(fileID,'1');
fclose(fileID);
end

