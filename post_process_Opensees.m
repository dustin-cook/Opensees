% Run Truss Tcl file in opensees
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = 'test';

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

%% Load element force data
for i = 1:length(element.id)
    ele_force(i,:) = max(dlmread([output_dir filesep ['element_force_' num2str(element.id(i)) '.txt']],' '));
end

% Add element forces to element table
element.Pmax = max([ele_force(:,1),ele_force(:,4)],[],2);
element.Vmax = max([ele_force(:,2),ele_force(:,5)],[],2);
element.Mmax = max([ele_force(:,3),ele_force(:,6)],[],2);

% % Omit elements that are rigid
% filter = (element.ele_id == 1) | (element.ele_id == 2);
% element(filter,:) = [];

%% Caclulate Element Capacity
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_id,:);
    [ ele_temp(i,:) ] = fn_element_capacity( ele, ele_prop );
end
element = ele_temp;

%% Calculate the DCR
element.DCR_P = element.Pmax ./ element.Pn;
element.DCR_V = element.Vmax ./ element.Vn_aci;
element.DCR_M = element.Mmax ./ element.Mn_aci;
DCR_max = max([element.DCR_P; element.DCR_V; element.DCR_M]);

% Save element table
writetable(element,[output_dir filesep 'element.csv'])

% Plot DCR view
s = element.node_1;
t = element.node_2;
G = graph(s,t);
for i = 1:max([element.node_1;element.node_2])
    if sum(node.id == i) == 0
        x(i) = 0;
        y(i) = 0;
    else
        x(i) = node.x(node.id == i);
        y(i) = node.y(node.id == i);
    end
end
% moment
s_break = element.node_1(element.DCR_M >= 1);
t_break = element.node_2(element.DCR_M >= 1);
H = plot(G,'XData',x,'YData',y);
highlight(H,s_break,t_break,'EdgeColor','red')
xlabel('Base (ft)')
xlabel('Height (ft)')
plot_name = 'DCR_view_moment';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

% Shear
s_break = element.node_1(element.DCR_V >= 1);
t_break = element.node_2(element.DCR_V >= 1);
H = plot(G,'XData',x,'YData',y);
highlight(H,s_break,t_break,'EdgeColor','red')
xlabel('Base (ft)')
xlabel('Height (ft)')
plot_name = 'DCR_view_shear';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

% Axial
s_break = element.node_1(element.DCR_P >= 1);
t_break = element.node_2(element.DCR_P >= 1);
H = plot(G,'XData',x,'YData',y);
highlight(H,s_break,t_break,'EdgeColor','red')
xlabel('Base (ft)')
xlabel('Height (ft)')
plot_name = 'DCR_view_axial';
fn_format_and_save_plot( plot_dir, plot_name, 4 )

%% Load Period data
periods = dlmread([output_dir filesep 'period.txt']);
T.x = periods(1);
T.z = periods(2);

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
        [ c1.(dirs_ran(i)), c2.(dirs_ran(i)) ] = fn_c_factors( site_class, num_stories, T.(dirs_ran(i)), haz_class, DCR_max );
    else
        c1.(dirs_ran(i)) = 1;
        c2.(dirs_ran(i)) = 1;
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
max_accel_all_nodes_x = max(abs(node.accel_x_abs),[],2)/386;
max_accel_profile_x = [max_accel_all_nodes_x(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    max_accel_profile_x(i+1) = max(max_accel_all_nodes_x(story.nodes_on_slab{i}));
end

% Z direction
if strcmp(model.dimension,'3D')
    max_accel_all_nodes_z = max(abs(node.accel_z_abs),[],2)/386;
    max_accel_profile_z = [max_accel_all_nodes_z(1), zeros(1,length(story.id))];
    for i = 1:length(story.id)
        max_accel_profile_z(i+1) = max(max_accel_all_nodes_z(story.nodes_on_slab{i}));
    end
end

%% Calculate Displacement Profile
% x direction
max_disp_all_nodes_x = max(abs(node.disp_x),[],2);
max_disp_profile_x = [max_disp_all_nodes_x(1), zeros(1,length(story.id))];
ave_disp_profile_x = [max_disp_all_nodes_x(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    max_disp_profile_x(i+1) = max(max_disp_all_nodes_x(story.nodes_on_slab{i}));
    ave_disp_profile_x(i+1) = mean(max_disp_all_nodes_x(story.nodes_on_slab{i}));
end

% Z direction
if strcmp(model.dimension,'3D')
    max_disp_all_nodes_z = max(abs(node.disp_z),[],2);
    max_disp_profile_z = [max_disp_all_nodes_z(1), zeros(1,length(story.id))];
    ave_disp_profile_z = [max_disp_all_nodes_z(1), zeros(1,length(story.id))];
    for i = 1:length(story.id)
        max_disp_profile_z(i+1) = max(max_disp_all_nodes_z(story.nodes_on_slab{i}));
        ave_disp_profile_z(i+1) = mean(max_disp_all_nodes_z(story.nodes_on_slab{i}));
    end
end

%% Torsional Amplification Check
TAR_x = max_disp_profile_x(2:end) ./ ave_disp_profile_x(2:end);
if sum(TAR_x > 1.2) > 0 && analysis.nonlinear == 0
    warning('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
end

if strcmp(model.dimension,'3D')
    TAR_z = max_disp_profile_z(2:end) ./ ave_disp_profile_z(2:end);
    if sum(TAR_z > 1.2) > 0 && analysis.nonlinear == 0
        warning('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
    end
end

% Amplify Displacements and Forces by Max TAR
if strcmp(model.dimension,'2D') && analysis.nonlinear == 0
    TAR_max = max(TAR_x);
    max_accel_profile_x = max_accel_profile_x*TAR_max;
    max_disp_profile_x = max_disp_profile_x*TAR_max;
end

%% Calculate Drift Profile
% X direction
max_drift_profile_x = zeros(length(story.id),1);
for i = 1:length(story.id)
    if i == 1
        nodal_drifts_x = node.disp_x(story.nodes_on_slab{i},:)/story.story_ht(i);
    else
        nodal_drifts_x = (node.disp_x(story.nodes_on_slab{i},:) - node.disp_x(story.nodes_on_slab{i-1},:))/story.story_ht(i);
    end
    max_drifts_all_slab_x = max(abs(nodal_drifts_x),[],2);
    max_drift_profile_x(i) = max(max_drifts_all_slab_x);
end

% Z direction
if strcmp(model.dimension,'3D')
    max_drift_profile_z = zeros(length(story.id),1);
    for i = 1:length(story.id)
        if i == 1
            nodal_drifts_z = node.disp_z(story.nodes_on_slab{i},:)/story.story_ht(i);
        else
            nodal_drifts_z = (node.disp_z(story.nodes_on_slab{i},:) - node.disp_z(story.nodes_on_slab{i-1},:))/story.story_ht(i);
        end
        max_drifts_all_slab_z = max(abs(nodal_drifts_z),[],2);
        max_drift_profile_z(i) = max(max_drifts_all_slab_z);
    end
end

%% Display Results
% Plot Roof Displacement
figure
hold on
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.disp_x(end,:),'DisplayName','X Direction')
if strcmp(model.dimension,'3D')
    plot((1:length(eq.z))*ground_motion.z.eq_dt,node.disp_z(end,:),'DisplayName','Z Direction')
end
xlabel('time (s)')
ylabel('Roof Displacemnet (in)')
plot_name = 'Roof_Displacemnet';
fn_format_and_save_plot( plot_dir, plot_name, 1)

% Plot Roof Acceleration
figure
hold on
xlabel('time (s)')
ylabel('Roof Acceleration (g)')
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.accel_x_abs(end,:)/386,'DisplayName','Roof')
plot((1:length(eq.x))*ground_motion.x.eq_dt,node.accel_x_abs(1,:)/386,'DisplayName','Ground')
plot_name = 'Roof_Acceleration_x';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

if strcmp(model.dimension,'3D')
    figure
    hold on
    xlabel('time (s)')
    ylabel('Roof Acceleration (g)')
    plot((1:length(eq.z))*ground_motion.z.eq_dt,node.accel_z_abs(end,:)/386,'DisplayName','Roof')
    plot((1:length(eq.z))*ground_motion.z.eq_dt,node.accel_z_abs(1,:)/386,'DisplayName','Ground')
    plot_name = 'Roof_Acceleration_z';
    fn_format_and_save_plot( plot_dir, plot_name, 1 )
end

% Plot Acceleration Profile
figure
hold on
plot(max_accel_profile_x,[0;story.id],'DisplayName','X Direction')
if strcmp(model.dimension,'3D')
    plot(max_accel_profile_z,[0;story.id],'DisplayName','Z Direction')
end
xlabel('PFA (g)')
ylabel('Story')
plot_name = 'Accel_Profile';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Displacement Profile
figure
hold on
plot(max_disp_profile_x,[0;story.id],'DisplayName','X Direction')
c0 = 1.3;
targ_disp = Sd.x*c1.x*c2.x*c0;
plot([0,targ_disp],[0,num_stories],'--r','DisplayName','Target Displacement')
if strcmp(model.dimension,'3D')
    plot(max_disp_profile_z,[0;story.id],'DisplayName','Z Direction')
end
xlabel('Displacement (in)')
ylabel('Story')
plot_name = 'Disp_Profile';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Drift Profile
figure
hold on
plot(max_drift_profile_x,story.id,'DisplayName','X Direction')
if strcmp(model.dimension,'3D')
    plot(max_drift_profile_z,story.id,'DisplayName','Z Direction')
end
xlabel('Story Drift')
ylabel('Story')
plot_name = 'Drift_Profile';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Save Data
save([output_dir filesep 'post_process_data'])
