function [ ] = fn_define_recorders( analysis_id )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Write Recorder File
file_name = ['Analysis' filesep analysis_id filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

% Define recorders
fprintf(fileID,'# Define recorders \n');
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_disp_x.txt -node 1 2 3 4 -dof 1 disp \n', analysis_id);
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_disp_y.txt -node 1 2 3 4 -dof 2 disp \n', analysis_id);
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_reaction_x.txt -node 1 2 3 4 -dof 1 reaction \n', analysis_id);
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_reaction_y.txt -node 1 2 3 4 -dof 2 reaction \n', analysis_id);
fprintf(fileID,'recorder Element -file Analysis/%s/element_1_force.txt -ele 1 force \n', analysis_id);
fprintf(fileID,'recorder Element -file Analysis/%s/element_2_force.txt -ele 2 force \n', analysis_id);
fprintf(fileID,'recorder Element -file Analysis/%s/element_3_force.txt -ele 3 force \n', analysis_id);

% Close File
fclose(fileID);

end

