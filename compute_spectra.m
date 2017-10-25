%% Method To calculate SDOP Response Spectra
clear
close
clc

%% Initailize Paramaters
% import hbr.str_eng.fn_sdof_th
dt = 0.01;
T = 0.05:0.05:3;
damp_ratio = 0.0;

%% Load in EQ data
load(['ground_motions' filesep 'data_ground_accel_th.mat']);
% load(['data_ground_motion' filesep 'data_atc_63_spectra.mat']);

%% Run SDOF response spectra\
eq_list = fieldnames(pga_th_g);
num_eq = numel(eq_list);
for i = 1:1
    ag = pga_th_g.(eq_list{i});
    for j = 1:length(T)
        [psuedoAccelerationTH, ~, accelerationTH] = fn_sdof_th(T(j), damp_ratio, ag, dt);
        max_sa(i,j) = max(abs(accelerationTH));
        max_sa_psuedo(i,j) = max(abs(psuedoAccelerationTH));    
    end
end

% med_spectra = median(max_sa,1);
% med_spectra_psuedo = median(max_sa_psuedo,1);
% med_pre_calc_spectra = median(sa_full_table(5:end,2:end),2);

%% Create Mean Response Spectra
hold on
grid on
% for i = 2:length(sa_full_table(1,:))
%     plot(sa_full_table(5:end,1),sa_full_table(5:end,i),'b','lineWidth',0.5)
% end
for i = 1:length(max_sa_psuedo(:,1))
    plot(T,max_sa(i,:),'b','lineWidth',0.5);
    plot(T,max_sa_psuedo(i,:),'k','lineWidth',0.5);
end
% plot(sa_full_table(5:end,1),med_pre_calc_spectra,'r','lineWidth',2.5)
% plot(T,mean_spectra,'c','lineWidth',2.5)
% plot(T,med_spectra_psuedo,'k','lineWidth',2)
