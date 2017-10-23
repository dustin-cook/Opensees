function [ ] = fn_define_recorders( analysis_id, nodes, ele )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% nodal_str = strjoin(num2str(nodes'),',');

% Write Recorder File
file_name = ['Analysis' filesep analysis_id filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

% Define recorders
fprintf(fileID,'# Define recorders \n');
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_disp_x.txt -node %s -dof 1 disp \n', analysis_id, num2str(nodes'));
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_disp_y.txt -node %s -dof 2 disp \n', analysis_id, num2str(nodes'));
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_reaction_x.txt -node %s -dof 1 reaction \n', analysis_id, num2str(nodes'));
fprintf(fileID,'recorder Node -file Analysis/%s/nodal_reaction_y.txt -node %s -dof 2 reaction \n', analysis_id, num2str(nodes'));
for i=1:length(ele)
    fprintf(fileID,'recorder Element -file Analysis/%s/element_%d_force.txt -ele %d force \n', analysis_id, ele(i), ele(i));
end

% Movie Recorders
fprintf(fileID,'# display displacement shape of the column \n');
fprintf(fileID,'recorder display "Displaced shape" 10 10 500 500 -wipe \n');
fprintf(fileID,'prp 200. 50. 1; \n');
fprintf(fileID,'vup 0 1 0; \n');
fprintf(fileID,'vpn 0 0 1; \n');
fprintf(fileID,'display 1 5 40 \n');

% Close File
fclose(fileID);

end

