clear all
close all
clc

%% Notes
% Roof Nodes: 6170, 6171, 6176, 6177

%% Inputs
dirs_ran = 'x';
% raw_accel_location = 'outputs\ICBS_model_3D_fixed\NDP_1\opensees_data\nodal_accel_6170.xml';
raw_accel_location = 'outputs\ICBS_model_3D_fixed\NDP_1\asce_41_data\node_TH_6170.mat';
% gm_location = 'ground_motions\ICSB_1979\chan_11_accel_50.tcl';
eq_dt = 0.01;
dt_scale = 0.01;

% %% Grab raw accels from opensees analysis and send to spectra tool
% [ node_accel_raw ] = fn_xml_read(raw_accel_location);
% node_accel_raw = node_accel_raw';
% 
% % Scale acceleration EQ (linear interpolation) based on uniform time step points
% eq = load(gm_location);
% eq_length = length(eq);
% eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
% eq_timespace_new = linspace(eq_dt,eq_length*eq_dt,eq_length*(eq_dt/dt_scale));
% eq_analysis_timespace = node_accel_raw(1,:);
% if strcmp(dirs_ran,'x')
%     node_accel_interp = interp1(eq_analysis_timespace,node_accel_raw(2,:),eq_timespace_new);
% else
%    node_accel_interp = interp1(eq_analysis_timespace,node_accel_raw(3,:),eq_timespace_new);
% end
% eq_interp = interp1(eq_timespace,eq',eq_timespace_new);
% 
% % scale to g and make absolute
% abs_accel_2_save = node_accel_interp/386 + eq_interp;

% filtered data
load(raw_accel_location)
abs_accel_2_save = nd_TH.(['accel_' dirs_ran '_abs_TH']);

% save 2 spectra tool inputs
save_dir = 'C:\Users\DustinCook\Desktop\Repos\toolbelt\Spectra Tool\inputs\ICSB_analysis';
fileID = fopen([save_dir filesep 'roof_accel_' dirs_ran '.txt.'],'w');
fprintf(fileID,'%s \n',dt_scale);
for i = 1:length(abs_accel_2_save)
    fprintf(fileID,'%s \n',abs_accel_2_save(i));
end
fclose(fileID);