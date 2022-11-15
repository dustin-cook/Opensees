function [ ] = fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import ida.fn_run_gm_ida
import ida.fn_run_gm_ida_sa_stripe

% Clear existing data or look through outputs to see if GM have already been run
gms2run = gm_set_table;
outputs_dir = [main_dir '/' 'IDA' ];
if analysis.clear_existing_data && isdir(outputs_dir)
    rmdir(outputs_dir, 's')
elseif exist(outputs_dir,'dir')
    if analysis.run_sa_stripes
        files = dir([outputs_dir filesep 'Summary Data' filesep 'GM_*']);
        final_sa = ['Sa_' strrep(num2str(analysis.sa_stripes(end)),'.','_')];
        for f = 1:length(files)
            if exist([outputs_dir filesep 'Summary Data' filesep files(f).name filesep final_sa filesep 'summary_results.mat'],'file')
                set_id = str2double(regexp(files(f).name,'(?<=_)\d+(?=_)','match'));
                pair = str2double(files(f).name(end));
                gms2run(gms2run.set_id == set_id & gms2run.pair == pair,:) = [];
            end
        end
    else
        files = dir([outputs_dir filesep 'GM_*']);
        for f = 1:length(files)
            if exist([outputs_dir filesep files(f).name filesep 'gm_complete.txt'],'file')
    %             completed_scales = dir([outputs_dir filesep files(f).name filesep 'Scale_*']);
    %             for s = 1:length(completed_scales)
    %                 rmdir([outputs_dir filesep files(f).name filesep completed_scales(s).name], 's')
    %             end
                set_id = str2double(regexp(files(f).name,'(?<=_)\d+(?=_)','match'));
                pair = str2double(files(f).name(end));
                gms2run(gms2run.set_id == set_id & gms2run.pair == pair,:) = [];
            end
        end
    end
end

% Set up Parallel Workers
if analysis.run_parallel
    parpool; 
end

% Loop through each ground motion
tim_start = tic;
if analysis.run_sa_stripes
    if analysis.run_parallel
        parfor gms = 1:height(gms2run)
            % Loop through each scale of the GM
            fn_run_gm_ida_sa_stripe(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gms, ida_results, tcl_dir, main_dir)
        end
    else
        for gms = 1:height(gms2run)
            % Loop through each scale of the GM
            fn_run_gm_ida_sa_stripe(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gms, ida_results, tcl_dir, main_dir)
        end
    end
else
    if analysis.run_parallel
        parfor gms = 1:height(gms2run)
            % Loop through each scale of the GM
            fn_run_gm_ida(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gms, ida_results, tcl_dir, main_dir)
        end
    else
        for gms = 1:height(gms2run)
            % Loop through each scale of the GM
            fn_run_gm_ida(analysis, model, story, element, node, hinge, joint, gm_set_table, gms2run, gms, ida_results, tcl_dir, main_dir)
        end
    end
end
tim_elapsed = toc(tim_start);
fprintf('IDA finished with a run time of %4.2f seconds \n', tim_elapsed)

delete(gcp('nocreate')) % End Any Parallel Process

end

