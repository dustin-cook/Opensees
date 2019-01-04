function [ ] = main_build_model( model, analysis, output_dir )
% Description: Main script facilitating building the model databases for
% the Opensees analysis

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
import build_model.fn_build_mdof
import build_model.fn_build_sdof

%% Begin Method
% Select Model Type
if analysis.model_type == 2 % MDOF Model
    fn_build_mdof( model, analysis, output_dir )
elseif analysis.model_type == 1 % SDOF Model
    fn_build_sdof( model, output_dir )
else
    error('Invalid Model Type')
end

end

