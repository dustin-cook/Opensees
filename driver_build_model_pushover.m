%% Clear the Workspace
clear
close
clc
fclose('all');

%% Description: Method to build an Opensees model and run a ASCE 41-17 teir 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 11/14/2022
% Inputs:

% Outputs:

% Assumptions:

%% User Inputs
analysis.id = 'Pushover_Model_LP'; % Unique String ID of the analysis for it to create its own directory

% Basic Analysis Options
analysis.nonlinear = 1; % Nonlinear
analysis.nonlinear_type = 'lumped'; % lumped or fiber
analysis.eq_lat_load_factor = 1;
analysis.type = 2; % Pushover
analysis.pushover_direction = 'x'; % {'x', '-x', 'z', '-z'} % Direction the pushover will run
analysis.dead_load = 1; % Dead Load Factor
analysis.live_load = 0.2; % Live Load Factor
analysis.accidental_torsion = 0; % add accidental torsion

% Model Options
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.stories_nonlinear_low = 0; % all stories at or below this story to be elastic (0 = all nonlinear)
analysis.elastic_beams = 0; % 0 = beams can be nonlinear (default), 1 = beams are assumed to be elastic
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)
analysis.fiber_walls = 0; % 0 = Lumped plasticity walls, 1 = Fiber walls, 2 = MVLEM
analysis.hinge_stiff_mod = 10; % Scale up stiffnes of hinges for a lumped plasticiy model. n value from Ibarra paper.
analysis.joint_model = 1; % 0 = centerline (almost), 1 = ASCE 41 implicit model, 2 = joint 3D, 3 = rigid (via beam/columns)
analysis.joint_explicit = 0; % 0 = rigid, 1 = model joint nonlinearity (could automate this based on first assessment of joints)
analysis.additional_elements = 1; % this is the leaning column

% Damping Options
analysis.damp_ratio = 0.03; % Critical damping ratio
analysis.damping = 'rayleigh'; % rayleigh, modal, or simple

% Recorder Options
analysis.simple_recorders = 0;
analysis.suppress_outputs = 1;
analysis.play_movie = 0;
analysis.hinge_group_length = 10;
analysis.write_xml = 1; % Write and read opensees out files as xml files (0 = .txt files, which is currently broken, on purpose)

% Other Analysis Options
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.opensees_SP = 0; % 0 = Standard OpenSees; 1 = OpenseesSP
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call

%% Initial Setup
% Import packages
import build_model.main_build_model
import opensees.write_tcl.*
import opensees.main_eigen_analysis

% Pull Model Data
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
num_models = height(model_data);

%% Initiate Analysis
for m = 1:num_models % run for each model    
    %% Load basic model data
    model = model_data(m,:);
    analysis.model_id = model.id;
    fprintf('Running Model %i of %i: %s\n', m, num_models, model.name{1})

    % Pull in database of available models
    model_table = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);

    % Select Model for the analysis 
    model = model_table(model_table.id == analysis.model_id,:);

    % Create Analysis Directory
    analysis.out_dir = ['outputs' filesep model.name{1} filesep analysis.id];
    mkdir(analysis.out_dir)

    %% Build CSV representation of model
    % Define CSV Directories
    analysis.model_dir = [analysis.out_dir filesep 'model_data'];
    
    % Connect elements to nodes and define joints and hinges
    main_build_model( model, analysis, [] )

    % Load Model Data from CSV files
    node = readtable([analysis.model_dir filesep 'node.csv'],'ReadVariableNames',true);
    element = readtable([analysis.model_dir filesep 'element.csv'],'ReadVariableNames',true);
    story = readtable([analysis.model_dir filesep 'story.csv'],'ReadVariableNames',true);
    joint = readtable([analysis.model_dir filesep 'joint.csv'],'ReadVariableNames',true);
    hinge = readtable([analysis.model_dir filesep 'hinge.csv'],'ReadVariableNames',true);

    %% Write TCL files
    % Create TCL Directory
    % TCL or Opensees does not like filesep command on windows, therefore must manually define forward slash seperators
    write_dir_opensees = [strrep(analysis.out_dir,'\','/') '/opensees_data']; 
    mkdir(write_dir_opensees)

    % Factor loads based on user input
    element.gravity_load =   analysis.dead_load*element.dead_load   + analysis.live_load*element.live_load;
    element.gravity_load_1 = analysis.dead_load*element.dead_load_1 + analysis.live_load*element.live_load_1;
    element.gravity_load_2 = analysis.dead_load*element.dead_load_2 + analysis.live_load*element.live_load_2;

    % Define Hinge Group (for recorders)
    hinge_grouping = [];
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        hinge.group = zeros(height(hinge),1);
        hinge_grouping = 1:analysis.hinge_group_length:(height(hinge)+1);
        if hinge_grouping(end) < (height(hinge)+1)
            hinge_grouping = [hinge_grouping, height(hinge)+1];
        end
        for hg = 1:(length(hinge_grouping)-1)
            group_range = hinge.id >= hinge_grouping(hg) & hinge.id < hinge_grouping(hg+1);
            hinge.group(group_range) = hg;
        end
    end

    % Write TCL model file
    [ ~ ] = fn_define_model( write_dir_opensees, node, element, joint, hinge, analysis, model.dimension, story, [], model );
    
    % Run Eigen Analysis to calc the first mode period
    [ model ] = main_eigen_analysis( model, analysis );
    
    % Write TCL files for recorders and loads
    fn_define_recorders( write_dir_opensees, model.dimension, node, element, joint, hinge, analysis );
    fn_define_loads( write_dir_opensees, analysis, node, model.dimension, story, element, joint, [], model);
end
