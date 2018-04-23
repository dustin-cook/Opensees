function [ ] = fn_define_recorders( output_dir, analysis, nodes )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% nodal_str = strjoin(num2str(nodes'),',');

% Write Recorder File
file_name = [output_dir, filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

% Define recorders
fprintf(fileID,'recorder Node -file %s/nodal_disp_x.txt -node %s -dof 1 disp \n', output_dir, num2str(nodes));
fprintf(fileID,'recorder Node -file %s/nodal_disp_y.txt -node %s -dof 2 disp \n', output_dir, num2str(nodes));
fprintf(fileID,'recorder Node -file %s/nodal_reaction_x.txt -node %s -dof 1 reaction \n', output_dir, num2str(nodes));
fprintf(fileID,'recorder Node -file %s/nodal_reaction_y.txt -node %s -dof 2 reaction \n', output_dir, num2str(nodes));
fprintf(fileID,'recorder Node -file %s/nodal_accel_x.txt -node %s -dof 1 accel \n', output_dir, num2str(nodes));
fprintf(fileID,'recorder Node -file %s/nodal_accel_y.txt -node %s -dof 2 accel \n', output_dir, num2str(nodes));
if strcmp(analysis.dims,'3D')
    fprintf(fileID,'recorder Node -file %s/nodal_disp_z.txt -node %s -dof 3 disp \n', output_dir, num2str(nodes));
    fprintf(fileID,'recorder Node -file %s/nodal_reaction_z.txt -node %s -dof 3 reaction \n', output_dir, num2str(nodes));
    fprintf(fileID,'recorder Node -file %s/nodal_accel_z.txt -node %s -dof 3 accel \n', output_dir, num2str(nodes));
end

% for i=1:length(ele)
%     fprintf(fileID,'recorder Element -file %s/element_%d_force.txt -ele %d force \n', output_dir, ele(i), ele(i));
% end

if analysis.type == 3
    % Movie Recorders
    fprintf(fileID,'recorder display "Displaced shape" 10 10 500 500 -wipe \n');
    fprintf(fileID,'prp 200. 50. 50.; \n');
    fprintf(fileID,'vup 0 1 0; \n');
    if strcmp(analysis.dims,'2D')
        fprintf(fileID,'vpn 0 0 1; \n');
    else
        fprintf(fileID,'vpn 0.4 0.25 1; \n');
    end
%     fprintf(fileID,'viewWindow -1000 1000 -1000 1000 \n');
    fprintf(fileID,'display 1 5 40 \n');
end

% Close File
fclose(fileID);

end

