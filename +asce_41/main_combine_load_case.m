function [ ] = main_combine_load_case( analysis )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Load Analysis Data
for i = 1:length(analysis.case_list)
    read_dir = [analysis.out_dir filesep analysis.case_list{i}];
    load([read_dir filesep 'story_analysis.mat']);
    load([read_dir filesep 'element_analysis.mat']);
    story_cases.(analysis.case_list{i}) = story;
    element_cases.(analysis.case_list{i}) = element;
end

%% Calculate Envelopes
% DCR
element.DCR_raw_max_M = max([element_cases.load_case_1.DCR_raw_max_M,element_cases.load_case_2.DCR_raw_max_M],[],2);
element.DCR_raw_max_V = max([element_cases.load_case_1.DCR_raw_max_V,element_cases.load_case_2.DCR_raw_max_V],[],2);
element.DCR_raw_max_P = max([element_cases.load_case_1.DCR_raw_max_P,element_cases.load_case_2.DCR_raw_max_P],[],2);
element.DCR_raw_max_all = max([element_cases.load_case_1.DCR_raw_max_all,element_cases.load_case_1.DCR_raw_max_all],[],2);
element.DCR_max_M = max([element_cases.load_case_1.DCR_max_M,element_cases.load_case_2.DCR_max_M],[],2);
element.DCR_max_V = max([element_cases.load_case_1.DCR_max_V,element_cases.load_case_2.DCR_max_V],[],2);
element.DCR_max_P = max([element_cases.load_case_1.DCR_max_P,element_cases.load_case_2.DCR_max_P],[],2);
element.DCR_max_all = max([element_cases.load_case_1.DCR_max_all,element_cases.load_case_1.DCR_max_all],[],2);

% EDP Profiles
story.max_accel_x = max([story_cases.load_case_1.max_accel_x,story_cases.load_case_2.max_accel_x],[],2);
story.max_disp_x = max([story_cases.load_case_1.max_disp_x,story_cases.load_case_2.max_disp_x],[],2);
story.max_drift_x = max([story_cases.load_case_1.max_drift_x,story_cases.load_case_2.max_drift_x],[],2);
if isfield(story,'max_drift_z')
    story.max_accel_z = max([story_cases.load_case_1.max_accel_z,story_cases.load_case_2.max_accel_z],[],2);
    story.max_disp_z = max([story_cases.load_case_1.max_disp_z,story_cases.load_case_2.max_disp_z],[],2);
    story.max_drift_z = max([story_cases.load_case_1.max_drift_z,story_cases.load_case_2.max_drift_z],[],2);
end

%% Save Data
write_dir = [analysis.out_dir filesep 'asce_41_data'];
save([write_dir filesep 'story_analysis.mat'],'story')
save([write_dir filesep 'element_analysis.mat'],'element')

end

