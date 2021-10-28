function [] = fn_run_gm_ida_sa_stripe(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gm_idx, ida_results, tcl_dir, main_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_main_IDA
import ida.fn_postprocess_ida
import ida.fn_postprocess_ida_gen

% Run for each scale factor
for i = 1:length(analysis.sa_stripes)
    % Define Ground Motion Data
    ground_motion.x = gms2run(gm_idx,:);
    ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
    ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
    if analysis.run_z_motion
        ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
        ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
        ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};
    end
    
    % Load spectral info and scale ground motion as geomean of pair
    spectra_table = readtable([ground_motion.x.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
    sa_gm_x = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1));
    if analysis.run_z_motion
        spectra_table = readtable([ground_motion.z.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
%         sa_gm_z = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(2));
        sa_gm_z = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1)); % use the first mode period for both directions to be consistent with USGS geomean
    else
        % Find the gm pair
        gm_alt = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
        gm_alt.eq_dir = {['ground_motions' '/' analysis.gm_set '/' gm_alt.eq_name{1}]};
        spectra_table = readtable([gm_alt.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
        sa_gm_z = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1));
    end
    if strcmp(analysis.scale_method,'geomean')
        sa_gm_geomean = sqrt(sa_gm_x*sa_gm_z);
        scale_factor = analysis.sa_stripes(i) / sa_gm_geomean;
    elseif strcmp(analysis.scale_method,'maxdir')
        sa_gm_maxdir = max(sa_gm_x,sa_gm_z);
        scale_factor = analysis.sa_stripes(i) / sa_gm_maxdir;
    end
    ground_motion.x.sa = sa_gm_x*scale_factor;
    if analysis.run_z_motion
        ground_motion.z.sa = sa_gm_z*scale_factor;
    end

    % Define Scale Factor
    analysis.ground_motion_scale_factor = scale_factor;
    
    % Create Output Directories
    ida_opensees_dir = [main_dir '/' 'IDA' '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair) '/' 'Sa_' strrep(num2str(analysis.sa_stripes(i)),'.','_')];
    if ~exist(ida_opensees_dir,'dir')
        mkdir(ida_opensees_dir)
    end

    ida_summary_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair) '/' 'Sa_' strrep(num2str(analysis.sa_stripes(i)),'.','_')];
    if ~exist(ida_summary_dir,'dir')
        mkdir(ida_summary_dir)
    end

    % Run Opensees
    if analysis.run_ida
        fprintf('Running GM ID: %i Pair: %i for an Sa of %4.2f g\n', gms2run.set_id(gm_idx), gms2run.pair(gm_idx), analysis.sa_stripes(i))
        tim_start = tic;
        [exit_status] = fn_main_IDA(analysis, model, story, element, node, hinge, joint, ground_motion, tcl_dir, ida_opensees_dir, ida_summary_dir);
        tim_elapsed = toc(tim_start);
        if exit_status == 1
            fprintf('Failed GM ID: %i Pair: %i for an Sa of %4.2f g\n', gms2run.set_id(gm_idx), gms2run.pair(gm_idx), analysis.sa_stripes(i))
        else
            fprintf('Successful run time of %4.2f seconds for GM ID: %i Pair: %i for an Sa of %4.2f g\n', tim_elapsed, gms2run.set_id(gm_idx), gms2run.pair(gm_idx), analysis.sa_stripes(i))
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
            fn_postprocess_ida(analysis, model, story, element, node, hinge, joint, ground_motion, ida_opensees_dir, ida_summary_dir)
        end
    end
end

end

