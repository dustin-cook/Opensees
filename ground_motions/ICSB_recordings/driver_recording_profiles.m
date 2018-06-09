clear 
close all
clc

%% Pull in recoding data and create EDP profiles
for i = 1:13
    accel_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_accel' filesep 'chan_' num2str(i) '_accel.tcl']);
    disp_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_disp' filesep 'chan_' num2str(i) '_disp.tcl']);
end

edp.floor = [0;1;2;3;4;5;6];
edp.ew_max_accel = [max(accel_TH.chan_13);max(accel_TH.chan_6);inf;max(accel_TH.chan_5);inf;inf;max(accel_TH.chan_4)];
edp.ns_max_accel = [max([accel_TH.chan_11;accel_TH.chan_10]);max([accel_TH.chan_7;accel_TH.chan_8;accel_TH.chan_9]);inf;inf;inf;inf;max([accel_TH.chan_1;accel_TH.chan_2;accel_TH.chan_3])];
ew_max_disp_rel = [max(disp_TH.chan_13);max(disp_TH.chan_6);inf;max(disp_TH.chan_5);inf;inf;max(disp_TH.chan_4)]/2.54; % Convert from cm to in (should do this earlier)
ns_max_disp_rel = [max([disp_TH.chan_11;disp_TH.chan_10]);max([disp_TH.chan_7;disp_TH.chan_8;disp_TH.chan_9]);inf;inf;inf;inf;max([disp_TH.chan_1;disp_TH.chan_2;disp_TH.chan_3])]/2.54;
edp.ew_max_disp = ew_max_disp_rel - ew_max_disp_rel(1);
edp.ns_max_disp = ns_max_disp_rel - ns_max_disp_rel(1);

edp_table = struct2table(edp);

writetable(edp_table,[pwd filesep 'recorded_edp_profile.csv'])