function [ dirs_ran, ground_motion ] = fn_load_gm_data( ground_motion_seq, gm_table )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
end

