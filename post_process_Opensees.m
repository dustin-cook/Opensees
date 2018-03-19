% Run Truss Tcl file in opensees
clear
close
clc

%% Load Analysis and Model parameters
analysis.model_id = 4;
analysis.gm_id = 4;
analysis.name = 'ASCE_41_LRHA';

gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
ground_motion = gm_table(gm_table.id == analysis.gm_id,:);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs/' model.name{1} '/' analysis.name];

%% Load outputs
load([output_dir filesep 'analysis_data.mat'])
eq_x = load([ground_motion.eq_dir{1} filesep ground_motion.eq_name_x{1}]);
eq_y = load([ground_motion.eq_dir{1} filesep ground_motion.eq_name_y{1}]);
eq_z = load([ground_motion.eq_dir{1} filesep ground_motion.eq_name_z{1}]);

% Nodal Displacement (in)
node.disp_x = dlmread([output_dir filesep 'nodal_disp_x.txt'],' ')';
node.disp_y = dlmread([output_dir filesep 'nodal_disp_y.txt'],' ')';
node.disp_z = dlmread([output_dir filesep 'nodal_disp_z.txt'],' ')';

% % Nodal Reaction (k)
% node.reac_x = dlmread([output_dir filesep 'nodal_reaction_x.txt'],' ')';
% node.reac_y = dlmread([output_dir filesep 'nodal_reaction_y.txt'],' ')';
% node.reac_z = dlmread([output_dir filesep 'nodal_reaction_z.txt'],' ')';

% Nodal Acceleration (in/s2)
node.accel_x_rel = dlmread([output_dir filesep 'nodal_accel_x.txt'],' ')';
node.accel_y_rel = dlmread([output_dir filesep 'nodal_accel_y.txt'],' ')';
node.accel_z_rel = dlmread([output_dir filesep 'nodal_accel_z.txt'],' ')';
node.accel_x_abs = node.accel_x_rel + ones(length(node.id),1)*eq_x'*386;
node.accel_y_abs = node.accel_y_rel + ones(length(node.id),1)*eq_y'*386;
node.accel_z_abs = node.accel_z_rel + ones(length(node.id),1)*eq_z'*386;

% % Element Forces
% for k = 1:length(element.id)
%     ele_force = dlmread([output_dir filesep 'element_' num2str(k) '_force.txt'],' ');
%     element.fx1{k} = ele_force(:,1)';
%     element.fy1{k} = ele_force(:,2)';
%     element.fz1{k} = ele_force(:,3)';
%     element.mx1{k} = ele_force(:,4)';
%     element.my1{k} = ele_force(:,5)';
%     element.mz1{k} = ele_force(:,6)';
%     element.fx2{k} = ele_force(:,7)';
%     element.fy2{k} = ele_force(:,8)';
%     element.fz2{k} = ele_force(:,9)';
%     element.mx2{k} = ele_force(:,10)';
%     element.my2{k} = ele_force(:,11)';
%     element.mz2{k} = ele_force(:,12)';
% end

%% Display Results
plot_dir = [output_dir filesep 'plots'];

% Plot Roof Displacement
figure
hold on
plot((1:length(eq_x))*ground_motion.eq_dt,node.disp_x(end,:),'DisplayName','X Direction')
plot((1:length(eq_z))*ground_motion.eq_dt,node.disp_z(end,:),'DisplayName','Z Direction')
xlabel('time (s)')
ylabel('Roof Displacemnet (in)')
plot_name = 'Roof_Displacemnet.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1)

% Plot Roof Acceleration
figure
hold on
xlabel('time (s)')
ylabel('Roof Acceleration (g)')
plot((1:length(eq_x))*ground_motion.eq_dt,node.accel_x_abs(end,:)/386,'DisplayName','Roof')
plot((1:length(eq_x))*ground_motion.eq_dt,node.accel_x_abs(1,:)/386,'DisplayName','Ground')
plot_name = 'Roof_Acceleration_x.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

figure
hold on
xlabel('time (s)')
ylabel('Roof Acceleration (g)')
plot((1:length(eq_z))*ground_motion.eq_dt,node.accel_z_abs(end,:)/386,'DisplayName','Roof')
plot((1:length(eq_z))*ground_motion.eq_dt,node.accel_z_abs(1,:)/386,'DisplayName','Ground')
plot_name = 'Roof_Acceleration_z.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Acceleration Profile
max_accel_all_nodes_x = max(abs(node.accel_x_abs),[],2)/386;
max_accel_all_nodes_z = max(abs(node.accel_z_abs),[],2)/386;
max_accel_profile_x = [max_accel_all_nodes_x(1), zeros(1,length(story.id))];
max_accel_profile_z = [max_accel_all_nodes_x(1), zeros(1,length(story.id))];
for i = 1:length(story.id)
    max_accel_profile_x(i+1) = max(max_accel_all_nodes_x(story.nodes_on_slab{i}));
    max_accel_profile_z(i+1) = max(max_accel_all_nodes_z(story.nodes_on_slab{i}));
end
figure
hold on
plot(max_accel_profile_x,[0;story.id],'DisplayName','X Direction')
plot(max_accel_profile_z,[0;story.id],'DisplayName','Z Direction')
xlabel('PFA (g)')
ylabel('Story')
plot_name = 'Accel_Profile.fig';
fn_format_and_save_plot( plot_dir, plot_name, 1 )

% Plot Drift Profile
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



