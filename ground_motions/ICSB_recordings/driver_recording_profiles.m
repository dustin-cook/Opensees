clear 
close all
clc

%% Pull in recoding data and create EDP profiles
for i = 1:13
    accel_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_accel' filesep 'chan_' num2str(i) '_accel.tcl']);
    disp_TH.(['chan_' num2str(i)]) = dlmread([pwd filesep 'chan_' num2str(i) '_disp' filesep 'chan_' num2str(i) '_disp.tcl']);
end

% EDP TH
record_edp.accel_TH_ground.x = accel_TH.chan_13;
record_edp.accel_TH_second.x = accel_TH.chan_6;
record_edp.accel_TH_roof.x = accel_TH.chan_4;
record_edp.accel_TH_ground.z = accel_TH.chan_11;
record_edp.accel_TH_second.z = accel_TH.chan_8;
record_edp.accel_TH_second_east.z = accel_TH.chan_9;
record_edp.accel_TH_roof.z = accel_TH.chan_2;
record_edp.accel_TH_roof_east.z = accel_TH.chan_3;
record_edp.disp_TH_ground.x = (disp_TH.chan_13(1:5800)-disp_TH.chan_13(1:5800))/2.54;
record_edp.disp_TH_second.x = (disp_TH.chan_6(1:5800)-disp_TH.chan_13(1:5800))/2.54;
record_edp.disp_TH_roof.x = (disp_TH.chan_4(1:5800)-disp_TH.chan_13(1:5800))/2.54;
record_edp.disp_TH_ground.z = (disp_TH.chan_11(1:5800)-disp_TH.chan_11(1:5800))/2.54;
record_edp.disp_TH_second.z = (disp_TH.chan_8(1:5800)-disp_TH.chan_11(1:5800))/2.54;
record_edp.disp_TH_second_east.z = (disp_TH.chan_9(1:5800)-disp_TH.chan_11(1:5800))/2.54;
record_edp.disp_TH_roof.z = (disp_TH.chan_2(1:5800)-disp_TH.chan_11(1:5800))/2.54;
record_edp.disp_TH_roof_east.z = (disp_TH.chan_3(1:5800)-disp_TH.chan_11(1:5800))/2.54;

% Max Profiles
record_edp.max_accel.x = [max(abs(accel_TH.chan_13));max(abs(accel_TH.chan_6));nan;max(abs(accel_TH.chan_5));nan;nan;max(abs(accel_TH.chan_4))];
record_edp.max_accel.z = [max(abs(accel_TH.chan_11));max(abs([accel_TH.chan_7;accel_TH.chan_8;accel_TH.chan_9]));nan;nan;nan;nan;max(abs([accel_TH.chan_1;accel_TH.chan_2;accel_TH.chan_3]))];
record_edp.max_disp.x = [max(abs(disp_TH.chan_13-disp_TH.chan_13));max(abs(disp_TH.chan_6-disp_TH.chan_13(1:5806)));nan;max(abs(disp_TH.chan_5-disp_TH.chan_13(1:5808)));nan;nan;max(abs(disp_TH.chan_4-disp_TH.chan_13(1:5806)))]/2.54; % Convert from cm to in (should do this earlier)
record_edp.max_disp.z = [max(abs(disp_TH.chan_11-disp_TH.chan_11));max(abs([disp_TH.chan_7(1:5808)-disp_TH.chan_11;disp_TH.chan_8-disp_TH.chan_11(1:5807);disp_TH.chan_9(1:5808)-disp_TH.chan_11]));nan;nan;nan;nan;max(abs([disp_TH.chan_1-disp_TH.chan_11(1:5806);disp_TH.chan_2-disp_TH.chan_11;disp_TH.chan_3-disp_TH.chan_11(1:5805)]))]/2.54;
record_edp.max_accel_center.x = [max(abs(accel_TH.chan_13));max(abs(accel_TH.chan_6));nan;max(abs(accel_TH.chan_5));nan;nan;max(abs(accel_TH.chan_4))];
record_edp.max_accel_center.z = [max(abs(accel_TH.chan_11));max(abs(accel_TH.chan_8));nan;nan;nan;nan;max(abs(accel_TH.chan_2))];
record_edp.max_disp_center.x = [max(abs(disp_TH.chan_13-disp_TH.chan_13));max(abs(disp_TH.chan_6-disp_TH.chan_13(1:5806)));nan;max(abs(disp_TH.chan_5-disp_TH.chan_13(1:5808)));nan;nan;max(abs(disp_TH.chan_4-disp_TH.chan_13(1:5806)))]/2.54; % Convert from cm to in (should do this earlier)
record_edp.max_disp_center.z = [max(abs(disp_TH.chan_11-disp_TH.chan_11));max(abs(disp_TH.chan_8-disp_TH.chan_11(1:5807)));nan;nan;nan;nan;max(abs(disp_TH.chan_2-disp_TH.chan_11))]/2.54;

save([pwd filesep 'recorded_edp_profile.mat'],'record_edp')