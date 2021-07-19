clear all
close all
clc

analysis.gm_set = 'FEMA_far_field';

gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);


% outputs_dir = [main_dir '/' 'IDA' ];
% files = dir([outputs_dir filesep 'GM_*']);
% for f = 1:length(files)
%     if exist([outputs_dir filesep files(f).name filesep 'gm_complete.txt'],'file')
% %             completed_scales = dir([outputs_dir filesep files(f).name filesep 'Scale_*']);
% %             for s = 1:length(completed_scales)
% %                 rmdir([outputs_dir filesep files(f).name filesep completed_scales(s).name], 's')
% %             end
%         set_id = str2double(regexp(files(f).name,'(?<=_)\d+(?=_)','match'));
%         pair = str2double(files(f).name(end));
%         gm_set_table(gm_set_table.set_id == set_id & gm_set_table.pair == pair,:) = [];
%     end
% end

for gms = 1:height(gm_set_table)

    ground_motion.x = gms2run(gm_idx,:);
    ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
    ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
    ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
    ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
    ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};

    % Load spectral info, define scale factor, and save Sa
    spectra_table = readtable([ground_motion.x.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
    sa_gm_x = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1));
    spectra_table = readtable([ground_motion.z.eq_dir{1} filesep 'spectra.csv'],'ReadVariableNames',true);
    sa_gm_z = interp1(spectra_table.period,spectra_table.psa_5,ida_results.period(1)); % use the first mode period for both directions to be consistent with USGS geomean
    sa_gm_geomean = geomean([sa_gm_x,sa_gm_z]);

    scale_factor = analysis.sa_stripes(i) / sa_gm_geomean;
    ground_motion.x.sa = sa_gm_x*scale_factor;
    if analysis.run_z_motion
        ground_motion.z.sa = sa_gm_z*scale_factor;
    end

end