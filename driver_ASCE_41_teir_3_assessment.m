%% Clear the Workspace
clear
close
clc

%% Description: Method to build an Opensees model and run a ASCE 41-17 teir
% 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% User Inputs (Think about changing this to a file read and command line execution)
analysis.model_id = 4;
analysis.proceedure = 'NDP'; % LDP or NDP or test
analysis.name = 'test'; % can remove this in lue of the procedure
analysis.dead_load = 1.0;
analysis.live_load = 1.0;
analysis.accidental_torsion = 0;

analysis.run_opensees = 1;
analysis.summit_SP = 0; % Write tcl files to be run on summit using OpenseesSP

analysis.gm_seq_id = 6; % Maybe also make this part ot the defualts or model?
analysis.pushover_drift = 0.02; % Maybe also make this part ot the defualts or model?
analysis.pushover_num_steps = 100; % Maybe also make this part ot the defualts or model?
analysis.pushover_direction = 'x'; % Maybe also make this part ot the defualts or model?
analysis.foundation = 1; % 2 = piles, 1 = fixed, 0 = pined % Maybe also make this part ot the defualts or model?
analysis.primary_node_offset = 0; % need to make this part of the model data or rework this input

%% Initial Setup
import asce_41.main_ASCE_41
import asce_41.fn_analysis_options

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Initiate Analysis
tic
main_ASCE_41( analysis )
toc

