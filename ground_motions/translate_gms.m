%% Method To calculate SDOP Response Spectra
clear
close
clc


%% Load in EQ data
load('data_ground_accel_th.mat');

%% Run SDOF response spectra\
eq_list = fieldnames(pga_th_g);
num_eq = numel(eq_list);
for i = 1:num_eq
    ag = pga_th_g.(eq_list{i});
    file_name = ['ATC_63_far_feild' filesep eq_list{i} '.tcl'];
    fileID = fopen(file_name,'w');
    fprintf(fileID,'%d \n',ag);
    fclose(fileID);
end

