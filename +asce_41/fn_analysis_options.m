function [ analysis ] = fn_analysis_options( analysis )
% Description: Defaults secondary analysis options for ASCE 41 Assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs: Analysis Data Structure

% Outputs: Analysis Data Structure

% Assumptions:

%% Begin Method
analysis.stories_nonlinear = inf; % Default to all modeling all stories as nonlinear when doing NDP
analysis.rigid_diaphram = 1; % Default the model to assume rigid diaphrams (0 = non-rigid assuption)

end

