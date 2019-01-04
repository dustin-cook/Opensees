% Run Post Analysis Logic for ASCE 41-17
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 6;
analysis.gm_id = 6;
analysis.name = 'test';
analysis.nonlinear = 0;

%% Import Packages
import asce_41.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
m_table.col = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.beam = readtable(['+asce_41' filesep 'linear_beam_m.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
load([output_dir filesep 'element_TH.mat'])
load([output_dir filesep 'model_analysis.mat'])
load([output_dir filesep 'element_analysis.mat'])
load([output_dir filesep 'hinge_analysis.mat'])
load([output_dir filesep 'story_analysis.mat'])
load([output_dir filesep 'gm_data.mat'])

% Linear Procedures
if analysis.nonlinear == 0
%     %% Caclulate Element Capacity and M factors
    for i = 1:length(element.id)
%         disp(['element ' num2str(i), ' out of ', num2str(length(element.id))])
        ele = element(i,:);
        ele_id = ele.ele_id;
        ele_prop = ele_prop_table(ele_id,:);
        ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
        
        % Calculate M Factors
        [ ele_temp(i,:) ] = fn_m_factors( m_table, ele, ele_prop );
    end
    element = ele_temp;

    %% Calculate the DCR for the shear calc
    [ element, model.DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'cp' ); %% NEED TO ITERATE ON DCR's HERE and maybe capacity as well ...

    % Perform calcs For each direction
    for i = 1:length(dirs_ran)
        if ~strcmp(dirs_ran{i},'y') % Update way I am doing this directional thing
        %% Caclulate C1 and C2 Factors
        model.num_stories = length(story.id);
        [ model.(['c1_' dirs_ran{i}]), model.(['c2_' dirs_ran{i}]), model.(['cm_' dirs_ran{i}]) ] = fn_c_factors( model.site_class{1}, model.num_stories, model.(['T1_' dirs_ran{i}]), model.(['hazus_class_' dirs_ran{i}]), model.DCR_raw_max );
        end
    end
    
    %% Amplify Element Forces ????HOW DO I DO THIS FOR 2 DIRECTIONS????
    element.Pmax_ASCE = element.Pmax*model.c1_x*model.c2_x;
    element.Vmax_ASCE = element.Vmax*model.c1_x*model.c2_x;
    element.Mmax_ASCE = element.Mmax*model.c1_x*model.c2_x;
    
    %% Torsional Amplification Check
    if strcmp(model.dimension,'3D')% for 3D linear analysis
        [ TAR_x, A_s_x ] = fn_tar_check( story.max_disp_x, story.ave_disp_x );
%         [ TAR_z, A_s_z ] = fn_tar_check( story.max_disp_z, story.ave_disp_z );
        element.Pmax_ASCE = element.Pmax_ASCE*max(A_s_x); % UPDATE TO WORK FOR WALL DIRECTIONS, AND WORK PER FLOOR
        element.Vmax_ASCE = element.Vmax_ASCE*max(A_s_x);
        element.Mmax_ASCE = element.Mmax_ASCE*max(A_s_x);
        story.max_disp_x_ASCE = story.max_disp_x*model.c1_x*model.c2_x .* A_s_x';
%         story.max_disp_z_ASCE = story.max_disp_z*c1.z*c2.z .* A_s_z';
        story.max_drift_x_ASCE = story.max_drift_x*model.c1_x*model.c2_x .* A_s_x';
%         story.max_drift_z_ASCE = story.max_drift_z*c1.z*c2.z .* A_s_z';
    else % For 2D linear analysis
        % Amplify Displacements and Forces by Max TAR
        story.max_accel_x = story.max_accel_x*model.tar;
        story.max_disp_x = story.max_disp_x*model.tar;
        element.Pmax_ASCE = element.Pmax_ASCE*model.tar;
        element.Vmax_ASCE = element.Vmax_ASCE*model.tar;
        element.Mmax_ASCE = element.Mmax_ASCE*model.tar;
        story.max_disp_x_ASCE = story.max_disp_x*model.c1_x*model.c2_x*model.tar;
        story.max_drift_x_ASCE = story.max_drift_x*model.c1_x*model.c2_x*model.tar;
    end
    
    %% Calculate the DCR
    [ element, model.DCR_raw_max ] = fn_calc_dcr( element, element_TH, 'cp' ); % UPDATE 2 directions

% Nonlinear Procedures
else 
    [ hinge ] = fn_accept_hinge( hinge, output_dir );
    story.max_disp_x_ASCE = story.max_disp_x;
    story.max_drift_x_ASCE = story.max_drift_x;
    if strcmp(model.dimension,'3D')
%         story.max_disp_z_ASCE = story.max_disp_z;
%         story.max_drift_z_ASCE = story.max_drift_z;
    end
end

%% Calculate Beam Column Strength Ratios
% for i =1:length(joint.id)
%    beam1 = max([element.Mn_pos(element.node_2 == joint.x_neg(i)),element.Mn_neg(element.node_2 == joint.x_neg(i))]); % Maximum of beam pos and neg nominal bending strength
%    beam2 = max([element.Mn_pos(element.node_1 == joint.x_pos(i)),element.Mn_neg(element.node_1 == joint.x_pos(i))]); 
%    column1 = min([element.Mn_pos(element.node_2 == joint.y_neg(i)),element.Mn_neg(element.node_2 == joint.y_neg(i))]); % Minimum of column pos and negative nominal moment strength
%    column2 = min([element.Mn_pos(element.node_1 == joint.y_pos(i)),element.Mn_neg(element.node_1 == joint.y_pos(i))]); 
%    joint.beam_strength(i) = sum([beam1,beam2]);
%    joint.column_strength(i) = sum([column1,column2]);
%    joint.col_bm_ratio(i) = joint.column_strength(i)/joint.beam_strength(i);
% end

%% Save Data
save([output_dir filesep 'hinge_analysis.mat'],'hinge')
save([output_dir filesep 'story_analysis.mat'],'story')
save([output_dir filesep 'model_analysis.mat'],'model')


