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
analysis.model_id = 2;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'Pushover'; % LDP or NDP or test
analysis.id = '1'; % ID of the analysis for it to create its own directory
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% dynamic analysis inputs
analysis.gm_seq_id = 16; % Maybe also make this part ot the defaults or model?

% Archetpye Inputs
analysis.additional_elements = 1; % this is the leaning column
analysis.eq_lat_load_factor = 1;

%% Initial Setup
import asce_41.main_ASCE_41

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Initiate Analysis
tic
main_ASCE_41( analysis )
toc

