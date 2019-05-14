clear all 
close all
clc

analysis.model_id = 11;
analysis.proceedure = 'NDP';
analysis.id = 2;
analysis.summit = 1;
analysis.run_ida = 1;
analysis.post_process_ida = 0;
analysis.gm_set = 'FEMA_far_field';

hazard.curve.rp = [43, 72, 224, 475, 975, 2475, 4975];%[43, 72, 224, 475, 975, 2475, 4975];
hazard.curve.pga = [0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];%[0.224, 0.308, 0.502, 0.635, 0.766, 0.946, 1.082];
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);
IDA_scale_factors = hazard.curve.pga ./ gm_median_pga;

model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);


for i = 1:length(IDA_scale_factors)
    for gm_idx = 1:height(gm_set_table)
        
        scale_factor = IDA_scale_factors(i);
        
        % Defin gms for this run
        ground_motion.x = gm_set_table(gm_idx,:);
        ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
        ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
        ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
        ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
        ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};

        % Iteration Parameters
        analysis.ground_motion_scale_factor = scale_factor;


        opensees_outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Scale_' num2str(analysis.ground_motion_scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];
        ida_outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(analysis.ground_motion_scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];
        mkdir(ida_outputs_dir)
        
        if exist([opensees_outputs_dir filesep 'hinge_analysis.mat'],'file')
            load([opensees_outputs_dir filesep 'hinge_analysis.mat'])
            load([opensees_outputs_dir filesep 'summary_results.mat'])
            load([opensees_outputs_dir filesep 'story_analysis.mat'])

            save([ida_outputs_dir filesep 'hinge_analysis.mat'],'hinge')
            save([ida_outputs_dir filesep 'summary_results.mat'],'summary')
            save([ida_outputs_dir filesep 'story_analysis.mat'],'story')
        end
    
    end
end