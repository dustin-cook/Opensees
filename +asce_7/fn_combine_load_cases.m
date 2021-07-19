function [ ] = fn_combine_load_cases( analysis, load_case_id )
% Description: Main script that post process an ELFP runs and combine loads
% from various load combos

% Created By: Dustin Cook
% Date Created: 7/9/2021

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'opensees_data'];
write_dir = [analysis.out_dir filesep 'asce_7_data'];
if ~exist(write_dir,'dir')
    fn_make_directory( write_dir )
end

% Load Analysis Data
% load([read_dir filesep 'model_analysis.mat'])
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'element_analysis.mat'])

% % Hard coded var that should be passed in as inputs
% Cd = 5.5;
% Ie = 1;
% 
% %% Caluclate Story Drifts
% disp = Cd * [0; story.ave_disp_x] / Ie;
% story.drift = abs(disp(2:end) - disp(1:(end-1))) ./ story.story_ht;

%% Combine demands from various load cases
if load_case_id == 1 % this is the first load case, just save the data
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
else % combine with previously run load cases
    % load previously saved demands
    story_combo = load([write_dir filesep 'story_analysis.mat']);
    element_combo = load([write_dir filesep 'element_analysis.mat']);

    % Merge element demands from most recent load case
    element.Pmax = max(element_combo.element.Pmax,element.Pmax);
    element.Pmin = min(element_combo.element.Pmin,element.Pmin);
    element.V = max(element_combo.element.V,element.V);
    element.Mpos = max(element_combo.element.Mpos,element.Mpos);
    element.Mneg = min(element_combo.element.Mneg,element.Mneg);

    % Merge story drift demands
    story.ave_disp_x = max(story_combo.story.ave_disp_x,story.ave_disp_x);
    story.max_disp_x = max(story_combo.story.max_disp_x,story.max_disp_x);
    
    % Overwrite saved load case
    save([write_dir filesep 'story_analysis.mat'],'story')
    save([write_dir filesep 'element_analysis.mat'],'element')
end

end
