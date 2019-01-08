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
analysis.model_id = 4;
analysis.proceedure = 'NDP'; % LDP or NDP or test
analysis.accidental_torsion = 0;

analysis.gm_seq_id = 8; % Maybe also make this part ot the defualts or model?

%% Initial Setup
import asce_41.main_ASCE_41
import asce_41.fn_analysis_options

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Initiate Analysis
tic
main_ASCE_41( analysis )
toc

