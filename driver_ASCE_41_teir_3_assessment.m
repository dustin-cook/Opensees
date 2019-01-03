% Description: Method to build an Opensees model and run a ASCE 41-17 teir
% 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:


% Clear the Workspace
clear
close
clc

%% User Inputs (Think about changing this to a file read and command line execution)
analysis.model_id = 6;
analysis.nonlinear = 0;
analysis.dead_load = 1.0;
analysis.live_load = 1.0;
analysis.accidental_torsion = 0;
analysis.primary_node_offset = 0; % need to make this part of the model data or rework this input
analysis.foundation = 1; % 2 = piles, 1 = fixed, 0 = pined, 

%% Initial Setup
import asce_41.main_ASCE_41
import asce_41.fn_analysis_options

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Initiate Analysis
tic
main_ASCE_41( analysis )
toc

