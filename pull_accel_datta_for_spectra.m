% Script to pull acceleration time histories from analysis and save them as
% spectra tool inputs

clear all
close all
clc

%% Inputs
analysis_dir = 'outputs\ICBS_model_3D_fixed\NDP_6\opensees_data';
model_dir = 'outputs\ICBS_model_3D_fixed\NDP_6\model_data';
spectra_dir = '..\toolbelt\Spectra Tool\inputs\ICSB_analysis';
center_x = 821;
center_z = 450;
dt = 0.01;

%% Load data
node = readtable([model_dir filesep 'node.csv'],'ReadVariableNames',true);

%% Find center nod
roof_nodes = node(node.y == max(node.y),:);
dist = sqrt((roof_nodes.x-center_x).^2 + (roof_nodes.z-center_z).^2);
[~,idx] = min(abs(dist));
node_center_roof = roof_nodes.id(idx);

%% Load accel TH data
load([analysis_dir filesep 'node_TH_' num2str(node_center_roof)])

%% Write to txt file
mkdir(spectra_dir)
% x dir
fileID = fopen([spectra_dir filesep 'roof_accel_x.txt'],'w');
fprintf(fileID,'%f \n', dt);
for i = 1:length(nd_TH.accel_x_abs_TH)
fprintf(fileID,'%f \n',nd_TH.accel_x_abs_TH(i));
end
fclose(fileID);
% z dir
fileID = fopen([spectra_dir filesep 'roof_accel_z.txt'],'w');
fprintf(fileID,'%f \n', dt);
for i = 1:length(nd_TH.accel_z_abs_TH)
fprintf(fileID,'%f \n',nd_TH.accel_z_abs_TH(i));
end
fclose(fileID);