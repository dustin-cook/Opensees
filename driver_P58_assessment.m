%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
% analysis.model_id = 4;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'P58'; 
analysis.id = '1'; % ID of the analysis for it to create its own directory
analysis.gm_set = 'FEMA_far_field';
analysis.run_z_motion = 0;

% Analysis options
analysis.summit = 0;
analysis.run_parallel = 0;
analysis.run_ida = 1;
analysis.post_process_ida = 1;
% analysis.create_fragilities = 0;
% analysis.plot_ida = 0;
% analysis.detialed_post_process = 0;
analysis.run_sa_stripes = 1;
% analysis.scale_increment = 0.25;
% analysis.sa_stripes = [0.2 0.4];
analysis.collapse_drift = 0.10; % Drift ratio at which we are calling the model essentially collapsed
analysis.clear_existing_data = 0;
analysis.general_ida = 0;
analysis.write_xml = 1;

% Nonlinear options
analysis.nonlinear = 1;
analysis.stories_nonlinear = inf;
analysis.stories_nonlinear_low = 0;
analysis.elastic_beams = 0;
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;

% Secondary options
analysis.dead_load = 1; % Dead load factor
analysis.live_load = 1; % live load factor
analysis.opensees_SP = 0;
analysis.type = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.03;
% analysis.run_eigen = 0;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 0;
analysis.play_movie = 0;
analysis.movie_scale = 10;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';
analysis.joint_explicit = 0;
analysis.simple_recorders = 1;
analysis.joint_model = 0;
analysis.additional_elements = 1;

% Site inputs
% Curt's thesis site
% site.lat = 33.996;
% site.lng = -118.162;
% Generic LA site
site.lat = 34.05;
site.lng = -118.25;

%% Import Packages
import ida.*
import asce_7.*
import build_model.main_build_model
import opensees.write_tcl.fn_define_model
import opensees.main_eigen_analysis
import usgs.*

model_ids = 1:24;
for m = 1:length(model_ids) % run for each model
    analysis.model_id = m;
%% Initial Setup
% Load basic model data
model_table = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

% Define read and write directories
main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
model_dir = [main_dir '/' 'model_data'];
if ~exist(model_dir,'dir')
    mkdir(model_dir)
end
tcl_dir = [main_dir '/' 'opensees_data'];
if ~exist(tcl_dir,'dir')
    mkdir(tcl_dir)
end
ELFP_model_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'model_data'];
% tcl_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'eigen_analysis'];
% asce41_dir = [main_dir '/' 'asce_41_data'];

% Load ground motion data
gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
gm_median_pga = median(gm_set_table.pga);

%% Build Model
analysis.out_dir = main_dir;
disp('Building Model ...')

% Create model tables
main_build_model( model, analysis, [] )

% Load Model Data
node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
element = readtable([model_dir filesep 'element.csv'],'readVariableNames',true);
joint = readtable([model_dir filesep 'joint.csv'],'readVariableNames',true);

% Write model tcl file
[ ~ ] = fn_define_model( tcl_dir, node, element, joint, hinge, analysis, model.dimension, story, [], model );

%% Run Eigen Analysis
% analysis.nonlinear = 0;
[ model ] = main_eigen_analysis( model, analysis );

%% Modify Model Data
% Set period variable
ida_results.period = model.T1_x;

% Factor Loads
[ element ] = fn_factor_loads( analysis, element, [] );

%% Pull Site Hazard Data from USGS
% [design_values] = fn_call_USGS_design_value_API(1, 'asce7-16', site.lat, site.lng, 'II', model.site_class);
afe = 1/475;
[sa_475] = fn_call_USGS_hazard_API('E2014', site.lat, site.lng, model.T1_x, model.vs30, afe);
analysis.sa_stripes = sa_475; % just run motions at the 475 year return period

%% Run Opensees Models
if analysis.run_ida || analysis.post_process_ida
    fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
end

%% Post Process Details of Specific GMs
% Collect EDP data
disp('Collecting EDP data for FEMA P-58 Assessment ...')
id = 0;
idr = [];
pfa = [];
for gm = 1:height(gm_set_table)
    gm_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm))];
    sa_folders = dir([gm_dir filesep 'Sa_*']);
    for s = 1:length(sa_folders)
        % Load data
        outputs_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm)) '/' sa_folders(s).name];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        story_file = [outputs_dir filesep 'story_analysis.mat'];
        if exist(outputs_file,'file') && exist(story_file,'file')
            load(outputs_file)
            load(story_file)
            
            id = id + 1;
            
            % Formulate the IDR table for the P-58 assessment - X direction
            idr.sa(id,1) = analysis.sa_stripes(s);
            idr.direction(id,1) = 1;
            idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
            for n = 1:height(story)
                idr.(['story_' num2str(n)])(id,1) = story.max_drift_x(n);
            end
            
            % Formulate the PFA table for the P-58 assessment - X direction
            pfa.sa(id,1) = analysis.sa_stripes(s);
            pfa.direction(id,1) = 1;
            pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
            pfa.floor_1(id,1) = summary.pga_x;
            for n = 1:height(story)
                pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_x(n);
            end
            
            id = id + 1;
            
            % Formulate the IDR table for the P-58 assessment - Z direction 
            if analysis.run_z_motion
                idr.sa(id,1) = analysis.sa_stripes(s);
                idr.direction(id,1) = 2;
                idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                for n = 1:height(story)
                    idr.(['story_' num2str(n)])(id,1) = story.max_drift_z(n);
                end
            else
                idr.sa(id,1) = analysis.sa_stripes(s);
                idr.direction(id,1) = 2;
                idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                for n = 1:height(story)
                    idr.(['story_' num2str(n)])(id,1) = story.max_drift_x(n);
                end
            end
            
            % Formulate the PFA table for the P-58 assessment - Z direction 
            if analysis.run_z_motion
                pfa.sa(id,1) = analysis.sa_stripes(s);
                pfa.direction(id,1) = 2;
                pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                pfa.floor_1(id,1) = summary.pga_z;
                for n = 1:height(story)
                    pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_z(n);
                end
            else
                pfa.sa(id,1) = analysis.sa_stripes(s);
                pfa.direction(id,1) = 2;
                pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                pfa.floor_1(id,1) = summary.pga_x;
                for n = 1:height(story)
                    pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_x(n);
                end
            end
            
            % Check for Collapse
            if summary.collapse > 0
                error('NEED TO EXPLICITLY HANDLE COLLAPSE CASES')
            end
        else
            error('NEED TO EXPLICITLY HANDLE MISSING STRIPES CASES')
        end
    end
end

% Convert to tables and sort
idr_table = struct2table(idr);
pfa_table = struct2table(pfa);
idr_table_sort = [];
pfa_table_sort = [];
for i = 1:length(sa_folders)
    for d = 1:2
        filt = idr_table.sa == analysis.sa_stripes(i) & idr_table.direction == d;
        idr_table_sect = idr_table(filt,:);
        idr_table_sort = [idr_table_sort; idr_table_sect];
        pfa_table_sect = idr_table(filt,:);
        pfa_table_sort = [pfa_table_sort; pfa_table_sect];
    end
end

% Write tables to save directory
write_dir = [main_dir '/' 'IDA' '/' 'EDPs'];
mkdir(write_dir)
writetable(idr_table_sort,[write_dir filesep 'idr.csv'])
writetable(pfa_table_sort,[write_dir filesep 'pfa.csv'])
end

