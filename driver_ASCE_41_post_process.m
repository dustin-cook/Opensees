% Run Post Analysis Logic for ASCE 41-17
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '11DL11LL';

%% Import Packages
import asce_41.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
m_table.col = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.beam = readtable(['+asce_41' filesep 'linear_beam_m.csv'],'ReadVariableNames',true);
load([output_dir filesep 'post_process_data.mat'])

%% Caclulate Element Capacity and M factors
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    [ ele ] = fn_element_capacity( ele, ele_prop );
    [ ele_temp(i,:) ] = fn_m_factors( m_table, ele, ele_prop );
end
element = ele_temp;


%% Calculate the DCR for the C1 calc
[ ~, DCR_max ] = fn_calc_dcr( element, 'cp', 1, 1, 'high' );

% Perform calcs For each direction
for i = 1:length(dirs_ran)
    %% Load Spectra and Calculate Sa
    spectra_table.(dirs_ran(i)) = readtable([ground_motion.(dirs_ran(i)).eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.(dirs_ran(i)).eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
    Sa.(dirs_ran(i)) = interp1(spectra_table.(dirs_ran(i)).period,spectra_table.(dirs_ran(i)).psa_5,model.(['T1_' dirs_ran(i)]));
    Sd.(dirs_ran(i)) = Sa.(dirs_ran(i))*386*(model.(['T1_' dirs_ran(i)])/(2*pi))^2;
    [ seismicity.(dirs_ran(i)) ] = fn_level_of_seismicity( spectra_table.(dirs_ran(i)) );

    %% Caclulate C1 and C2 Factors
    if analysis.nonlinear == 0
        [ c1.(dirs_ran(i)), c2.(dirs_ran(i)) ] = fn_c_factors( model.site_class{1}, length(story.id), model.(['T1_' dirs_ran(i)]), model.(['hazus_class_' num2str(i)]), DCR_max );
    else
        c1.(dirs_ran(i)) = 1;
        c2.(dirs_ran(i)) = 1;
    end
    
    %% Amplify Element Forces ????HOW DO I DO THIS FOR 2 DIRECTIONS????
    element.Pmax_ASCE = element.Pmax*c1.(dirs_ran(i))*c2.(dirs_ran(i));
    element.Vmax_ASCE = element.Vmax*c1.(dirs_ran(i))*c2.(dirs_ran(i));
    element.Mmax_ASCE = element.Mmax*c1.(dirs_ran(i))*c2.(dirs_ran(i));
end

%% Torsional Amplification Check
if strcmp(model.dimension,'3D') && analysis.nonlinear == 0 % for 3D linear analysis
    [ TAR_x ] = fn_tar_check( max_disp_profile_x, ave_disp_profile_x );
    [ TAR_z ] = fn_tar_check( max_disp_profile_z, ave_disp_profile_z );
elseif analysis.nonlinear == 0 % For 2D linear analysis
    % Amplify Displacements and Forces by Max TAR
    story.max_accel_x = story.max_accel_x*model.tar;
    story.max_disp_x = story.max_disp_x*model.tar;
    element.Pmax_ASCE = element.Pmax_ASCE*model.tar;
    element.Vmax_ASCE = element.Vmax_ASCE*model.tar;
    element.Mmax_ASCE = element.Mmax_ASCE*model.tar;
end

%% Calculate the DCR
[ element, ~ ] = fn_calc_dcr( element, 'cp', c1.x, c2.x, seismicity.x ); % UPDATE FOR SEISMICITY AND 2 directions

%% Save Data
% remove extra data
clear analysis

% save data
save([output_dir filesep 'ASCE_data'])


