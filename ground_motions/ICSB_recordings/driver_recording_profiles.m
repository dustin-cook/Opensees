clear 
close all
clc

%% Pull in recoding data and create EDP profiles
for i = 1:13
    accel_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_accel' filesep 'chan_' num2str(i) '_accel.tcl']);
    disp_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_disp' filesep 'chan_' num2str(i) '_disp.tcl']);
end

record_edp.max_accel.x = [max(accel_TH.chan_13);max(accel_TH.chan_6);nan;max(accel_TH.chan_5);nan;nan;max(accel_TH.chan_4)];
record_edp.max_accel.z = [max([accel_TH.chan_11;accel_TH.chan_10]);max([accel_TH.chan_7;accel_TH.chan_8;accel_TH.chan_9]);nan;nan;nan;nan;max([accel_TH.chan_1;accel_TH.chan_2;accel_TH.chan_3])];
ew_max_disp_rel = [max(disp_TH.chan_13);max(disp_TH.chan_6);nan;max(disp_TH.chan_5);nan;nan;max(disp_TH.chan_4)]/2.54; % Convert from cm to in (should do this earlier)
ns_max_disp_rel = [max([disp_TH.chan_11;disp_TH.chan_10]);max([disp_TH.chan_7;disp_TH.chan_8;disp_TH.chan_9]);nan;nan;nan;nan;max([disp_TH.chan_1;disp_TH.chan_2;disp_TH.chan_3])]/2.54;
record_edp.max_disp.x = ew_max_disp_rel - ew_max_disp_rel(1);
record_edp.max_disp.z = ns_max_disp_rel - ns_max_disp_rel(1);

save([pwd filesep 'recorded_edp_profile.mat'],'record_edp')