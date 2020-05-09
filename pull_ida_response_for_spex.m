%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
analysis.model_id = 11;
analysis.proceedure = 'NDP';
analysis.id = 'IDA_new';
analysis.gm_set = 'FEMA_far_field';
analysis.run_z_motion = 1;

%% Initial Setup
% Load basic model data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Define directories
ida_output_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id filesep 'IDA' filesep 'Summary Data'];

% Define hazard points to use
hazard_sa = [0.22, 0.328, 0.422, 0.492];

%% Loop through each ground motion and scale factors folder and grab peak response data
gm_dirs = dir([ida_output_dir filesep 'GM_*']);
id = 0;
for g = 1:length(gm_dirs)
    gm_dir = gm_dirs(g).name;
    scale_dirs = dir([ida_output_dir filesep gm_dir filesep 'Scale_*']);
    for s = 1:length(scale_dirs)
        scale_dir = scale_dirs(s).name;
        load([ida_output_dir filesep gm_dir filesep scale_dir filesep 'story_analysis.mat'])
        load([ida_output_dir filesep gm_dir filesep scale_dir filesep 'summary_results.mat'])
        
        % Pull response values and save in table
        id = id + 1;
        response.gm_id(id,1) = g;
        response.sa_x(id,1) = summary.sa_x;
        response.sa_z(id,1) = summary.sa_z;
        response.idr_x{id,1} = story.max_drift_x;
        response.idr_z{id,1} = story.max_drift_z;
        response.ridr_x(id,1) = max([story.residual_disp_x(1); story.residual_disp_x(2:end) - story.residual_disp_x(1:(end-1))] ./ story.story_ht);
        response.ridr_z(id,1) = max([story.residual_disp_z(1); story.residual_disp_z(2:end) - story.residual_disp_z(1:(end-1))] ./ story.story_ht);
    end
    
    % Find peak sa of each ground motion
    gm_ids(g) = g;
    peak_sa_x(g) = max(response.sa_x(response.gm_id == g));
    peak_sa_z(g) = max(response.sa_z(response.gm_id == g));
end

response_table = struct2table(response);

%% remove ground motions that do not get up tp hazard levels
gms_ids_to_remove = gm_ids(peak_sa_x < max(hazard_sa));
response_table(ismember(response_table.gm_id,gms_ids_to_remove),:) = [];
gms_to_use = unique(response_table.gm_id);

%% Go through each intensity and ground motion and interpolate for response at specified hazard points
for h = 1:length(hazard_sa)
    for g = 1:length(gms_to_use)
        gm = response_table(response_table.gm_id == gms_to_use(g),:);
        % Develop idr arrays
        idr_x = [];
        idr_z = [];
        for i = 1:height(gm)
            idr_x = [idr_x, gm.idr_x{i}];
            idr_z = [idr_z, gm.idr_z{i}];
        end
        interp_idr_x(g,:) = interp1(gm.sa_x,idr_x',hazard_sa(h));
        interp_idr_z(g,:) = interp1(gm.sa_x,idr_z',hazard_sa(h));
        interp_ridr_x(g) = interp1(gm.sa_x,gm.ridr_x,hazard_sa(h));
        interp_ridr_z(g) = interp1(gm.sa_z,gm.ridr_z,hazard_sa(h));
        
%         % Write SPEx Input strings for all gms inputs
%         fprintf('data.responses.per_intensity{%i}.all_gms{%i}.idr = [ %s ; %s ]; \n', h, g, num2str(interp_idr_x), num2str(interp_idr_z));
%         fprintf('data.responses.per_intensity{%i}.all_gms{%i}.ridr_max = [ %s ; %s ]; \n', h, g, num2str(interp_ridr_x), num2str(interp_ridr_z));
    end
    % Write SPEx Input strings for median inputs
    fprintf('data.responses.per_intensity{%i}.med{1}.idr = [ %s ]; \n', h, num2str(median(interp_idr_x)));
    fprintf('data.responses.per_intensity{%i}.med{2}.idr = [ %s ]; \n', h, num2str(median(interp_idr_z)));
    fprintf('data.responses.per_intensity{%i}.med{1}.ridr_max = [ %s ]; \n', h, num2str(median(interp_ridr_x)));
    fprintf('data.responses.per_intensity{%i}.med{2}.ridr_max = [ %s ]; \n', h, num2str(median(interp_ridr_z)));
end

