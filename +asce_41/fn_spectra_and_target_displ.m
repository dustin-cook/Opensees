function [ target_disp_in ] = fn_spectra_and_target_displ( model, story, ground_motion, direction )
% Description: Main script facilitating an ASCE 41 teir 3 assessment. 

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Import Packages
import asce_41.fn_cm
import asce_41.fn_target_disp

%% Begin Method
% Load Spectra and Calculate Sa
spectra_table = readtable([ground_motion.(direction).eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.(direction).eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
Sa = interp1(spectra_table.period,spectra_table.psa_5,model.(['T1_' direction]));

% Calculate Target Displacement
strength_ratio = model.DCR_raw_max; % We have no Vy for linear analysis, therefor use DCR max as a proxy for strength ratio
[ c_m ] = fn_cm( model.num_stories, model.(['T1_' direction]), model.(['hazus_class_' direction]) );
[ target_disp_in ] = fn_target_disp( strength_ratio, model.site_class{1}, story.(['mode_shape_' direction]), model.num_stories, model.(['T1_' direction]), Sa, c_m );

end

