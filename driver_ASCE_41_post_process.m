% Run Post Analysis Logic for ASCE 41-17
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 4;
analysis.gm_id = 8;
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
load([output_dir filesep 'element_TH.mat'])
load([output_dir filesep 'element_analysis.mat'])

% Linear Procedures
if analysis.nonlinear == 0
    %% Caclulate Element Capacity and M factors
    for i = 1:length(element.id)
        disp(['element ' num2str(i), ' out of ', num2str(length(element.id))])
        ele = element(i,:);
        ele_id = ele.ele_id;
        ele_prop = ele_prop_table(ele_id,:);
        ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
        
        % Calculate Element Capacity
        [ ele, element_TH.(['ele_' num2str(element.id(i))]) ] = fn_element_capacity( ele, ele_prop, ele_TH );
        
        % Caculate required development length and make sure there is enough
        [ ele.pass_aci_dev_length ] = fn_development_check( ele, ele_prop );
        
        % Calculate M Factors
        [ ele_temp(i,:) ] = fn_m_factors( m_table, ele, ele_prop );
    end
    element = ele_temp;

    %% Calculate the DCR for the shear calc
    [ element, DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'cp' );

    %% Caclulate Element Capacity and M factors
    for i = 1:length(element.id)
        ele = element(i,:);
        ele_id = ele.ele_id;
        ele_prop = ele_prop_table(ele_id,:);
        ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
        
        % Calculate Element Capacity
        [ ele, element_TH.(['ele_' num2str(element.id(i))]) ] = fn_element_capacity( ele, ele_prop, ele_TH );
        
        % Caculate required development length and make sure there is enough
        [ ele.pass_aci_dev_length ] = fn_development_check( ele, ele_prop );
        
        % Calculate M Factors
        [ ele_temp2(i,:) ] = fn_m_factors( m_table, ele, ele_prop );
    end
    element = ele_temp2;

    %% Calculate the DCR for the C1 calc
    [ element, DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'cp' );

    % Perform calcs For each direction
    for i = 1:length(dirs_ran)
        if ~strcmp(dirs_ran(i),'y') % Update way I am doing this directional thing
        %% Load Spectra and Calculate Sa
        spectra_table.(dirs_ran(i)) = readtable([ground_motion.(dirs_ran(i)).eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.(dirs_ran(i)).eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
        Sa.(dirs_ran(i)) = interp1(spectra_table.(dirs_ran(i)).period,spectra_table.(dirs_ran(i)).psa_5,model.(['T1_' dirs_ran(i)]));
        Sd.(dirs_ran(i)) = Sa.(dirs_ran(i))*386*(model.(['T1_' dirs_ran(i)])/(2*pi))^2;
        [ seismicity.(dirs_ran(i)) ] = fn_level_of_seismicity( spectra_table.(dirs_ran(i)) );

        %% Caclulate C1 and C2 Factors
        [ c1.(dirs_ran(i)), c2.(dirs_ran(i)) ] = fn_c_factors( model.site_class{1}, length(story.id), model.(['T1_' dirs_ran(i)]), model.(['hazus_class_' dirs_ran(i)]), DCR_raw_max );
        end
    end
    
    %% Amplify Element Forces ????HOW DO I DO THIS FOR 2 DIRECTIONS????
    element.Pmax_ASCE = element.Pmax*c1.x*c2.x;
    element.Vmax_ASCE = element.Vmax*c1.x*c2.x;
    element.Mmax_ASCE = element.Mmax*c1.x*c2.x;
    
    %% Torsional Amplification Check
    if strcmp(model.dimension,'3D')% for 3D linear analysis
        [ TAR_x, A_s_x ] = fn_tar_check( story.max_disp_x, story.ave_disp_x );
        [ TAR_z, A_s_z ] = fn_tar_check( story.max_disp_z, story.ave_disp_z );
        element.Pmax_ASCE = element.Pmax_ASCE*max(A_s_x); % UPDATE TO WORK FOR WALL DIRECTIONS, AND WORK PER FLOOR
        element.Vmax_ASCE = element.Vmax_ASCE*max(A_s_x);
        element.Mmax_ASCE = element.Mmax_ASCE*max(A_s_x);
        story.max_disp_x_ASCE = story.max_disp_x*c1.x*c2.x .* A_s_x';
        story.max_disp_z_ASCE = story.max_disp_z*c1.z*c2.z .* A_s_z';
        story.max_drift_x_ASCE = story.max_drift_x*c1.x*c2.x .* A_s_x';
        story.max_drift_z_ASCE = story.max_drift_z*c1.z*c2.z .* A_s_z';
    else % For 2D linear analysis
        % Amplify Displacements and Forces by Max TAR
        story.max_accel_x = story.max_accel_x*model.tar;
        story.max_disp_x = story.max_disp_x*model.tar;
        element.Pmax_ASCE = element.Pmax_ASCE*model.tar;
        element.Vmax_ASCE = element.Vmax_ASCE*model.tar;
        element.Mmax_ASCE = element.Mmax_ASCE*model.tar;
        story.max_disp_x_ASCE = story.max_disp_x*c1.x*c2.x*model.tar;
        story.max_drift_x_ASCE = story.max_drift_x*c1.x*c2.x*model.tar;
    end
    
    %% Calculate the DCR
    [ element, DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'cp' ); % UPDATE 2 directions

% Nonlinear Procedures
else 
    [ hinge ] = fn_accept_hinge( hinge, output_dir );
end

%% Save Data
% remove extra data
clear analysis

% save data
save([output_dir filesep 'ASCE_data'])


