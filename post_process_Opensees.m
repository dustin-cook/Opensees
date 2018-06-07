% Post Process Opensees Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '09DL';

%% Import Packages
import asce_41.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_id,:);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
plot_dir = [output_dir filesep 'plots'];

element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);

load([output_dir filesep 'analysis_data.mat'])

% Load ground motion data
dirs_ran = [];
if ground_motion_seq.eq_id_x~=0
    dirs_ran = [dirs_ran, 'x'];
    ground_motion.x = gm_table(gm_table.id == ground_motion_seq.eq_id_x,:);
end
if ground_motion_seq.eq_id_y~=0
    dirs_ran = [dirs_ran, 'y'];
    ground_motion.y = gm_table(gm_table.id == ground_motion_seq.eq_id_y,:);
end
if ground_motion_seq.eq_id_z~=0
    dirs_ran = [dirs_ran, 'z'];
    ground_motion.z = gm_table(gm_table.id == ground_motion_seq.eq_id_z,:);
end

% Load element force data
for i = 1:length(element.id)
    ele_force(i,:) = max(dlmread([output_dir filesep ['element_force_' num2str(element.id(i)) '.txt']],' '));
end

% Add element forces to element table
element.Pmax = max([ele_force(:,1),ele_force(:,4)],[],2);
element.Vmax = max([ele_force(:,2),ele_force(:,5)],[],2);
element.Mmax = max([ele_force(:,3),ele_force(:,6)],[],2);

% Load Period data
periods = dlmread([output_dir filesep 'period.txt']);
T.x = periods(1);
T.z = periods(2);

%% Caclulate Element Capacity
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    [ ele_temp(i,:) ] = fn_element_capacity( ele, ele_prop );
end
element = ele_temp;

%% Calculate the DCR
[ element, DCR_max ] = fn_calc_dcr( element );

% Perform calcs For each direction
for i = 1:length(dirs_ran)
    %% Load Spectra and Calculate Sa
    spectra_table.(dirs_ran(i)) = readtable([ground_motion.(dirs_ran(i)).eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.(dirs_ran(i)).eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
    Sa.(dirs_ran(i)) = interp1(spectra_table.(dirs_ran(i)).period,spectra_table.(dirs_ran(i)).psa_5,T.(dirs_ran(i)));
    Sd.(dirs_ran(i)) = Sa.(dirs_ran(i))*386*(T.(dirs_ran(i))/(2*pi))^2;

    %% Caclulate C1 and C2 Factors
    if analysis.nonlinear == 0
        site_class = model.site_class{1};
        num_stories = length(story.id);
        if i == 1
            haz_class = model.hazus_class_1;
        elseif i == 2
            haz_class = model.hazus_class_2;
        end
        DCR_max = 1;
        [ c1.(dirs_ran(i)), c2.(dirs_ran(i)) ] = fn_c_factors( site_class, num_stories, T.(dirs_ran(i)), haz_class, DCR_max );
    else
        c1.(dirs_ran(i)) = 1;
        c2.(dirs_ran(i)) = 1;
    end
    
    %% Amplify Element Forces
    if strcmp(model.dimension,'2D')
        element.Pmax_mod = element.Pmax*c1.x*c2.x;
        element.Vmax_mod = element.Vmax*c1.x*c2.x;
        element.Mmax_mod = element.Mmax*c1.x*c2.x;
    else
        warning('ADD LOGIC TO AMPLIFY ELEMENT FORCES BY C FACTORS FOR 3D')
    end

    %% Load and Read Outputs
    eq.(dirs_ran(i)) = load([ground_motion.(dirs_ran(i)).eq_dir{1} filesep ground_motion.(dirs_ran(i)).eq_name{1}]);

    % Nodal Displacement (in)
    node.(['disp_' dirs_ran(i)]) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*dlmread([output_dir filesep ['nodal_disp_' dirs_ran(i) '.txt']],' ')';

    % Nodal Acceleration (in/s2)
    node.(['accel_' dirs_ran(i) '_rel']) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*dlmread([output_dir filesep ['nodal_accel_' dirs_ran(i) '.txt']],' ')';
    node.(['accel_' dirs_ran(i) '_abs']) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*(node.(['accel_' dirs_ran(i) '_rel']) + ones(length(node.id),1)*eq.(dirs_ran(i))'*386);
end

%% Calculate Acceleration Profile
% x Direction
[ max_accel_profile_x ] = fn_calc_max_repsonse_profile( node.accel_x_abs, story, 0, 386 );

% Z direction
if strcmp(model.dimension,'3D') && strcmp(dirs_ran,'z')
    [ max_accel_profile_z ] = fn_calc_max_repsonse_profile( node.accel_z_abs, story, 0, 386 );
end

%% Calculate Displacement Profile
% x direction
[ max_disp_profile_x ] = fn_calc_max_repsonse_profile( node.disp_x, story, 0, 1 );
[ ave_disp_profile_x ] = fn_calc_max_repsonse_profile( node.disp_x, story, 1, 1 );

% Z direction
if strcmp(model.dimension,'3D') && strcmp(dirs_ran,'z')
    [ max_disp_profile_z ] = fn_calc_max_repsonse_profile( node.disp_z, story, 0, 1 );
    [ ave_disp_profile_z ] = fn_calc_max_repsonse_profile( node.disp_z, story, 1, 1 );
end

%% Torsional Amplification Check
if strcmp(model.dimension,'3D') && analysis.nonlinear == 0 % for 3D linear analysis
    [ TAR_x ] = fn_tar_check( max_disp_profile_x, ave_disp_profile_x );
    [ TAR_z ] = fn_tar_check( max_disp_profile_z, ave_disp_profile_z );
elseif analysis.nonlinear == 0 % For 2D linear analysis
    % Amplify Displacements and Forces by Max TAR
        TAR_max = model.tar;
        max_accel_profile_x = max_accel_profile_x*TAR_max;
        max_disp_profile_x = max_disp_profile_x*TAR_max;
        node.accel_x_rel = node.accel_x_rel*TAR_max;
        node.accel_x_abs = node.accel_x_abs*TAR_max;
        node.disp_x = node.disp_x*TAR_max;
        element.Pmax_mod = element.Pmax_mod*TAR_max;
        element.Vmax_mod = element.Vmax_mod*TAR_max;
        element.Mmax_mod = element.Mmax_mod*TAR_max;
end

%% Calculate Drift Profile
% X direction
[ max_drift_profile_x ] = fn_drift_profile( node.disp_x, story );

% Z direction
if strcmp(model.dimension,'3D')
    [ max_drift_profile_z ] = fn_drift_profile( node.disp_z, story );
end

%% Display Results
% Plot Roof Displacement
fn_plot_response_history( node.disp_x, eq.x, ground_motion.x.eq_dt, plot_dir, 'Roof Displacemnet X (in)' )

if strcmp(model.dimension,'3D')
    fn_plot_response_history( node.disp_z, eq.z, ground_motion.z.eq_dt, plot_dir, 'Roof Displacemnet Z (in)' )
end

% Plot Roof Acceleration
fn_plot_response_history( node.accel_x_abs/386, eq.x, ground_motion.x.eq_dt, plot_dir, 'Roof Acceleration X (g)' )

if strcmp(model.dimension,'3D')
    fn_plot_response_history( node.accel_z_abs/386, eq.z, ground_motion.z.eq_dt, plot_dir, 'Roof Acceleration Z (g)' )
end

% Plot Acceleration Profile
if strcmp(model.dimension,'3D')
    fn_plot_profile( max_accel_profile_x, [0;story.id], plot_dir, 'Acceleration Profile', 'PFA (g)', max_accel_profile_z )
    fn_plot_profile( max_disp_profile_x, [0;story.id], plot_dir, 'Displacement Profile', 'Displacement (in)', max_disp_profile_z )
    fn_plot_profile( max_drift_profile_x, story.id, plot_dir, 'Drift Profile', 'IDR', max_drift_profile_z )
else
    fn_plot_profile( max_accel_profile_x, [0;story.id], plot_dir, 'Acceleration Profile', 'PFA (g)' )
    fn_plot_profile( max_disp_profile_x, [0;story.id], plot_dir, 'Displacement Profile', 'Displacement (in)' )
    fn_plot_profile( max_drift_profile_x, story.id, plot_dir, 'Drift Profile', 'IDR' )
end

%% Calculate the DCR
element.DCR_P = element.Pmax_mod ./ element.Pn;
element.DCR_V = element.Vmax_mod ./ element.Vn_aci;
element.DCR_M = element.Mmax_mod ./ element.Mn_aci;

% Save element table
writetable(element,[output_dir filesep 'element.csv'])

%% Save Data
save([output_dir filesep 'post_process_data'])
