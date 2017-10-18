function [ ] = fn_define_loads( story_force_k, story_wieght_k, analysis_id )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Write Loads File
file_name = ['Analysis' filesep analysis_id filesep 'loads.tcl'];
fileID = fopen(file_name,'w');

% Define Gravity Loads (node id, axial, shear, moment)
fprintf(fileID,'# Define Gravity Loads (node id, axial, shear, moment) \n');
fprintf(fileID,'pattern Plain 1 Linear {  \n');
fprintf(fileID,'   load 3 0.0 %d 0.0 \n', story_wieght_k/2);
fprintf(fileID,'   load 4 0.0 %d 0.0 \n', story_wieght_k/2);
fprintf(fileID,'} \n');

% Define Load Pattern
fprintf(fileID,'# Define Load Pattern \n');
fprintf(fileID,'pattern Plain 2 Linear { \n');
fprintf(fileID,'  load 3 %d 0.0 0.0 \n', story_force_k);
fprintf(fileID,'} \n');

% Close File
fclose(fileID);

end

