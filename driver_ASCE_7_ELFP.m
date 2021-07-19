%% Clear the Workspace
clear
close
clc
fclose('all');

%% Description: Method to build an Opensees model and run a ASCE 41-17 teir 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% User Inputs (Think about changing this to a file read and command line execution)
analysis.model_id = 9;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'ELFP'; % LDP or NDP or test
analysis.id = '1'; % ID of the analysis for it to create its own directory

%% Initial Setup
import asce_7.main_ASCE_7
import build_model.main_build_model
import opensees.main_eigen_analysis

%% Secondary Fixed Inputs
analysis.run_opensees = 1;
analysis.run_eigen = 0;
analysis.opensees_SP = 0;
analysis.algorithm = 'KrylovNewton';
analysis.joint_model = 0;

analysis.suppress_outputs = 0;
analysis.simple_recorders = 0;
analysis.play_movie = 1;
analysis.movie_scale = 10;
analysis.summit = 0;
analysis.write_xml = 1;

%% Define Model and Run Eigen Analysis
model_table = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Create Analysis Directory
analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id];
fn_make_directory( analysis.out_dir )

% Build Model
disp('Building Model ...')
analysis.type = 4;
main_build_model( model, analysis, [] )

% Run Eigen Analysis
analysis.nonlinear = 0;
[ model ] = main_eigen_analysis( model, analysis );

% Write period to design sheet
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], model.T1_x, 'OpenseesOutput', 'C2')

%% Test Load Cases
% analysis.type_list =          [4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 = static loads
% analysis.nonlinear_list =     [0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
% analysis.drift_run_list =     [0]; % run exceptions for ELFP drift cases
% analysis.dead_load_list =     [1.2]; % Dead load factor for linear analysis
% analysis.live_out_load_list = [0.5]; % Live load factor on outer bays for linear analysis
% analysis.live_in_load_list =  [0.5]; % Live load factor on inner bays for linear analysis
% analysis.eq_lat_load_list =   [1]; % earthquake load factor for linear analysis
% analysis.eq_vert_load_list =  [0.2]; % earthquake vertical load factor for linear analysis
% 
% % Run the load case
% load_case_id = 'test';
% main_ASCE_7( analysis, load_case_id, model )

%% Drift Load Cases
analysis.type_list =          [4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 = static loads
analysis.nonlinear_list =     [0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
analysis.drift_run_list =     [1]; % run exceptions for ELFP drift cases
analysis.dead_load_list =     [0]; % Dead load factor for linear analysis
analysis.live_out_load_list = [0]; % Live load factor on outer bays for linear analysis
analysis.live_in_load_list =  [0]; % Live load factor on inner bays for linear analysis
analysis.eq_lat_load_list =   [1]; % earthquake load factor for linear analysis
analysis.eq_vert_load_list =  [0]; % earthquake vertical load factor for linear analysis

% Run the load case
tic
load_case_id = 'drift';
main_ASCE_7( analysis, load_case_id, model )
toc

% Write Displacement and story shear data to design sheet
asce_7_dir = [analysis.out_dir filesep load_case_id filesep 'asce_7_data'];
load([asce_7_dir filesep 'story_analysis.mat'])
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.max_disp_x', 'OpenseesOutput', 'C5')

%% Primary Force Load Cases
analysis.type_list =          [4, 4, 4, 4, 4, 4, 4, 4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 =  static loads
analysis.nonlinear_list =     [0, 0, 0, 0, 0, 0, 0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
analysis.drift_run_list =     [0, 0, 0, 0, 0, 0, 0, 0]; % run exceptions for ELFP drift cases
analysis.dead_load_list =     [1.4, 1.2, 1.2, 1.2, 1.2,  0.9, 1.2,  0.9]; % Dead load factor for linear analysis
analysis.live_out_load_list = [0,   1.6, 1.6, 0,   0.5,  0,   0.5,  0]; % Live load factor on outer bays for linear analysis
analysis.live_in_load_list =  [0,   1.6, 0,   1.6, 0.5,  0,   0.5,  0]; % Live load factor on outer bays for linear analysis
analysis.eq_lat_load_list =   [0,   0,   0,   0,   1,    1,  -1,   -1]; % earthquake load factor for linear analysis
analysis.eq_vert_load_list =  [0,   0,   0,   0,   0.2, -0.2, 0.2, -0.2]; % earthquake vertical load factor for linear analysis

% Run the load case
tic
load_case_id = 'forces';
main_ASCE_7( analysis, load_case_id, model )
toc

% Write Peak Bending Moments to Design Sheet
asce_7_dir = [analysis.out_dir filesep load_case_id filesep 'asce_7_data'];
load([asce_7_dir filesep 'element_analysis.mat'])
load([asce_7_dir filesep 'story_analysis.mat'])
bm_filt = strcmp(element.type,'beam');
col_filt = strcmp(element.type,'column');
inner_col = element.inner_bay == 1;
for s = 1:height(story)
    story_filt = element.story == s;
    story.beam_pos_bending(s) = max(element.Mpos(bm_filt & story_filt));
    story.beam_neg_bending(s) = max(element.Mneg(bm_filt & story_filt));
    story.ext_col_bending(s) = max(abs([element.Mpos(col_filt & story_filt & ~inner_col); element.Mneg(col_filt & story_filt & ~inner_col)]));
    story.int_col_bending(s) = max(abs([element.Mpos(col_filt & story_filt & inner_col); element.Mneg(col_filt & story_filt & inner_col)]));
end
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.max_disp_x', 'OpenseesOutput', 'C6')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.beam_pos_bending', 'OpenseesOutput', 'C7')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.beam_neg_bending', 'OpenseesOutput', 'C8')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.ext_col_bending', 'OpenseesOutput', 'C9')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.int_col_bending', 'OpenseesOutput', 'C10')

%% Overturning Load Case
analysis.type_list =          [4, 4, 4, 4, 4, 4]; % 1 = dynamic, 2 = pushover % 3 = static cyclic, 4 =  static loads
analysis.nonlinear_list =     [0, 0, 0, 0, 0, 0]; % 0 = linear, 1 = IMK Rotational Hinge, 2 = strain hardening hinges
analysis.drift_run_list =     [0, 0, 0, 0, 0, 0]; % run exceptions for ELFP drift cases
analysis.dead_load_list =     [1.4, 1.2, 1.2,  0.9,  1.2,  0.9]; % Dead load factor for linear analysis
analysis.live_out_load_list = [0,   1.6, 0.5,  0,    0.5,  0]; % Live load factor on outer bays for linear analysis
analysis.live_in_load_list =  [0,   1.6, 0.5,  0,    0.5,  0]; % Live load factor on inner bays for linear analysis
analysis.eq_lat_load_list =   [0,   0,   0.5,  0.5, -0.5, -0.5]; % earthquake load factor for linear analysis
analysis.eq_vert_load_list =  [0,   0,   0.2, -0.2,  0.2, -0.2]; % earthquake vertical load factor for linear analysis
    
% Run the load case
tic
load_case_id = 'overturning';
main_ASCE_7( analysis, load_case_id, model )
toc

% Write Peak Bending Moments to Design Sheet
asce_7_dir = [analysis.out_dir filesep load_case_id filesep 'asce_7_data'];
load([asce_7_dir filesep 'element_analysis.mat'])
load([asce_7_dir filesep 'story_analysis.mat'])
col_filt = strcmp(element.type,'column');
inner_col = element.inner_bay == 1;
for s = 1:height(story)
    story_filt = element.story == s;
    story.ext_col_max_axial(s) = max(element.Pmax(col_filt & story_filt & ~inner_col));
    story.ext_col_min_axial(s) = min(element.Pmin(col_filt & story_filt & ~inner_col));
    story.int_col_max_axial(s) = max(element.Pmax(col_filt & story_filt & inner_col));
    story.int_col_min_axial(s) = min(element.Pmin(col_filt & story_filt & inner_col));
end
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.ext_col_max_axial', 'OpenseesOutput', 'C11')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.ext_col_min_axial', 'OpenseesOutput', 'C12')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.int_col_max_axial', 'OpenseesOutput', 'C13')
xlswrite([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'], story.int_col_min_axial', 'OpenseesOutput', 'C14')
