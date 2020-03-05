function [ ] = fn_define_loads_ida( write_dir, ground_motion_scale_factor, dimension, ground_motion, g_unit )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import asce_7.*

% Write Loads File
file_name = [write_dir filesep 'loads.tcl'];
fileID = fopen(file_name,'w');

fprintf(fileID,'puts "Defining Loads ..." \n');


%% Dynamic Analysis
% Define Seismic Excitation Load
% timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>
% pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $vel0> <-fact $cFactor>
if isfield(ground_motion,'x')
    fprintf(fileID,'set dt %f \n',ground_motion.x.eq_dt);
%         fprintf(fileID,'puts "EQ X dt = $dt" \n');
    fprintf(fileID,'timeSeries Path 1 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.x.eq_dir{1}, ground_motion.x.eq_name{1}, g_unit);
    fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 -fact %f \n',ground_motion_scale_factor); 
end
if isfield(ground_motion,'z') && strcmp(dimension,'3D')
    fprintf(fileID,'set dt %f \n',ground_motion.z.eq_dt);
%         fprintf(fileID,'puts "EQ Z dt = $dt" \n');
    fprintf(fileID,'timeSeries Path 2 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.z.eq_dir{1}, ground_motion.z.eq_name{1}, g_unit);
    fprintf(fileID,'pattern UniformExcitation 4 3 -accel 2 -fact %f \n',ground_motion_scale_factor); 
end
if isfield(ground_motion,'y')
    fprintf(fileID,'set dt %f \n',ground_motion.y.eq_dt);
%         fprintf(fileID,'puts "EQ Y dt = $dt" \n');
    fprintf(fileID,'timeSeries Path 3 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.y.eq_dir{1}, ground_motion.y.eq_name{1}, g_unit);
    fprintf(fileID,'pattern UniformExcitation 5 2 -accel 3 -fact %f \n',ground_motion_scale_factor); 
end

fprintf(fileID,'puts "Define Load Complete" \n');

%% Close File
fclose(fileID);

end

