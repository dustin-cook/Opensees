function [ ] = main_build_model( model, analysis, ele_prop_table )
% Description: Main script facilitating building the model databases for
% the Opensees analysis

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
import build_model.fn_build_archetype
import build_model.fn_build_mdof
import build_model.fn_build_sdof

% Define Read Directory
read_dir = ['inputs' filesep 'models' filesep model.name{1}];

% Create Write Directory
write_dir = [analysis.out_dir filesep 'model_data'];
fn_make_directory( write_dir )

%% Begin Method
% Select Model Type
if analysis.model_type == 3 % Archetype Model
    fn_build_archetype( model, analysis, write_dir )
elseif analysis.model_type == 2 % MDOF Model
    fn_build_mdof( model, ele_prop_table, analysis, write_dir, read_dir )
elseif analysis.model_type == 1 % SDOF Model
    fn_build_sdof( model, analysis, write_dir )
else
    error('Invalid Model Type')
end

end

