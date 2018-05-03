% Run Truss Tcl file in opensees
clear
close
clc

%% Load Analysis and Model parameters
analysis.model_id = 1;
analysis.gm_id = 3;
analysis.name = 'rayleigh';

%% Load Analysis Data
gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_id,:);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs/' model.name{1} '/' analysis.name];
plot_dir = [output_dir filesep 'plots'];
load([output_dir filesep 'analysis_data.mat'])
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

%% Load Period data
periods = dlmread([output_dir filesep 'period.txt']);
T_1 = periods(1);
T_2 = periods(2);

%% Load Spectra and Calculate Sa
spectra_table_x = readtable([ground_motion.x.eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.x.eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
Sa_1 = interp1(spectra_table_x.period,spectra_table_x.psa_5,T_1);
if isfield(ground_motion,'z')
    spectra_table_z = readtable([ground_motion.z.eq_dir{1} filesep 'spectra_' erase(erase(ground_motion.z.eq_name{1},'.tcl'),'gm_') '.csv'],'ReadVariableNames',true);
    Sa_2 = interp1(spectra_table_z.period,spectra_table_z.psa_5,T_2);
end

%% Caclulate C1 and C2 Factors
if analysis.nonlinear == 0
    site_class = 'D';
    if strcmp(site_class,'A') || strcmp(site_class,'B')
        a = 130;
    elseif strcmp(site_class,'C')
        a = 90;
    else
        a = 60;
    end
    num_stories = length(story.id);
    if num_stories >= 3 || T_1 >= 1
        c_m = 1;
    else
        c_m = 0.9; % Need to update this to work in two directions and accept building types
    end
    DCR_max_1 = 4; % Placeholder for now until we get actual values
    DCR_max_2 = 4; % Placeholder for now until we get actual values
    u_strength_1 = max([DCR_max_1*c_m/1.5,1]);
    u_strength_2 = max([DCR_max_2*c_m/1.5,1]);
    c1.x = 1 + (u_strength_1-1)/(a*T_1^2);
    c2.x = 1 + (1/800)*((u_strength_1-1)/T_1)^2;
    c1.z = 1 + (u_strength_2-1)/(a*T_2^2);
    c2.z = 1 + (1/800)*((u_strength_2-1)/T_2)^2;
else
    c1.x = 1;
    c2.x = 1;
    c1.z = 1;
    c2.z = 1;
end
%% Load and Read Outputs
for i = 1:length(dirs_ran)
eq.(dirs_ran(i)) = load([ground_motion.(dirs_ran(i)).eq_dir{1} filesep ground_motion.(dirs_ran(i)).eq_name{1}]);

% Nodal Displacement (in)
node.(['disp_' dirs_ran(i)]) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*dlmread([output_dir filesep ['nodal_disp_' dirs_ran(i) '.txt']],' ')';

% Nodal Acceleration (in/s2)
node.(['accel_' dirs_ran(i) '_rel']) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*dlmread([output_dir filesep ['nodal_accel_' dirs_ran(i) '.txt']],' ')';
node.(['accel_' dirs_ran(i) '_abs']) = c1.(dirs_ran(i))*c2.(dirs_ran(i))*(node.(['accel_' dirs_ran(i) '_rel']) + ones(length(node.id),1)*eq.(dirs_ran(i))'*386);

end

%% EDP Profiles
% Acceleration
max_accel_all_nodes_x = max(abs(node.accel_x_abs),[],2)/386;
max_accel_all_nodes_z = max(abs(node.accel_z_abs),[],2)/386;
max_accel_profile_x = [max_accel_all_nodes_x(1), zeros(1,length(story.id))];
max_accel_profile_z = [max_accel_all_nodes_z(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    max_accel_profile_x(i+1) = max(max_accel_all_nodes_x(story.nodes_on_slab{i}));
    max_accel_profile_z(i+1) = max(max_accel_all_nodes_z(story.nodes_on_slab{i}));
end

% Displacement
max_disp_all_nodes_x = max(abs(node.disp_x),[],2);
max_disp_all_nodes_z = max(abs(node.disp_z),[],2);
max_disp_profile_x = [max_disp_all_nodes_x(1), zeros(1,length(story.id))];
max_disp_profile_z = [max_disp_all_nodes_z(1), zeros(1,length(story.id))];
ave_disp_profile_x = [max_disp_all_nodes_x(1), zeros(1,length(story.id))];
ave_disp_profile_z = [max_disp_all_nodes_z(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    max_disp_profile_x(i+1) = max(max_disp_all_nodes_x(story.nodes_on_slab{i}));
    max_disp_profile_z(i+1) = max(max_disp_all_nodes_z(story.nodes_on_slab{i}));
    ave_disp_profile_x(i+1) = mean(max_disp_all_nodes_x(story.nodes_on_slab{i}));
    ave_disp_profile_z(i+1) = mean(max_disp_all_nodes_z(story.nodes_on_slab{i}));
end

% Drift
max_drift_profile_x = zeros(length(story.id),1);
max_drift_profile_z = zeros(length(story.id),1);
for i = 1:length(story.id)
    if i == 1
        nodal_drifts_x = node.disp_x(story.nodes_on_slab{i},:)/story.story_ht(i);
        nodal_drifts_z = node.disp_z(story.nodes_on_slab{i},:)/story.story_ht(i);
    else
        nodal_drifts_x = (node.disp_x(story.nodes_on_slab{i},:) - node.disp_x(story.nodes_on_slab{i-1},:))/story.story_ht(i);
        nodal_drifts_z = (node.disp_z(story.nodes_on_slab{i},:) - node.disp_z(story.nodes_on_slab{i-1},:))/story.story_ht(i);
    end
    max_drifts_all_slab_x = max(abs(nodal_drifts_x),[],2);
    max_drifts_all_slab_z = max(abs(nodal_drifts_z),[],2);
    max_drift_profile_x(i) = max(max_drifts_all_slab_x);
    max_drift_profile_z(i) = max(max_drifts_all_slab_z);
end

%% Torsional Amplification Check
TAR_x = max_disp_profile_x(2:end) ./ ave_disp_profile_x(2:end);
TAR_z = max_disp_profile_z(2:end) ./ ave_disp_profile_z(2:end);

if sum(TAR_x > 1.2) > 0 || sum(TAR_z > 1.2) > 0 
    warning('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
end

%% Display Results
% Plot Roof Displacement
figure
hold on
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.disp_x(end,:),'DisplayName','X Direction')
plot((1:length(eq.z))*ground_motion.z.eq_dt,node.disp_z(end,:),'DisplayName','Z Direction')
xlabel('time (s)')
ylabel('Roof Displacemnet (in)')
plot_name = 'Roof_Displacemnet.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1)

% Plot Roof Acceleration
figure
hold on
xlabel('time (s)')
ylabel('Roof Acceleration (g)')
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.accel_x_abs(end,:)/386,'DisplayName','Roof')
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.accel_x_abs(1,:)/386,'DisplayName','Ground')
plot_name = 'Roof_Acceleration_x.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

figure
hold on
xlabel('time (s)')
ylabel('Roof Acceleration (g)')
plot((1:length(eq.z))*ground_motion.z.eq_dt,node.accel_z_abs(end,:)/386,'DisplayName','Roof')
plot((1:length(eq.z))*ground_motion.z.eq_dt,node.accel_z_abs(1,:)/386,'DisplayName','Ground')
plot_name = 'Roof_Acceleration_z.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Acceleration Profile
figure
hold on
plot(max_accel_profile_x,[0;story.id],'DisplayName','X Direction')
plot(max_accel_profile_z,[0;story.id],'DisplayName','Z Direction')
xlabel('PFA (g)')
ylabel('Story')
plot_name = 'Accel_Profile.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Drift Profile
figure
hold on
plot(max_drift_profile_x,story.id,'DisplayName','X Direction')
plot(max_drift_profile_z,story.id,'DisplayName','Z Direction')
xlabel('Story Drift')
ylabel('Story')
plot_name = 'Drift_Profile.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )


% Save Data
clear analysis
save([output_dir filesep 'analysis_data'])

