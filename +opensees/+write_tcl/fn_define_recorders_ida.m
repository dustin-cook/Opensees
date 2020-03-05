function [ ] = fn_define_recorders_ida( write_dir, dimension, node_ids, element_ids )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Define File Type to Write to
% if analysis.write_xml
    file_type = '-xml';
    file_ext = 'xml';
% else
%     file_type = '-file';
%     file_ext = 'txt';
% end

%     file_type = '-file';
%     file_ext = 'out';

% Write Recorder File
file_name = [write_dir, filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

fprintf(fileID,'puts "Defining Recorders ..." \n');

% Change the maximun number of files that can be open (only for windows
% machines, only matters for big models)
fprintf(fileID,'setMaxOpenFiles 2000 \n');


% Define Node recorders
for i = 1:length(node_ids)
    if strcmp(dimension,'2D')
            fprintf(fileID,'recorder Node %s %s/nodal_disp_%s.%s -time -node %i -dof 1 disp \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
            fprintf(fileID,'recorder Node %s %s/nodal_accel_%s.%s -time -node %i -dof 1 accel \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
            fprintf(fileID,'recorder Node %s %s/nodal_reaction_%s.%s -time -node %i -dof 1 3 reaction \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
    elseif strcmp(dimension,'3D')
            fprintf(fileID,'recorder Node %s %s/nodal_disp_%s.%s -time -node %i -dof 1 3 disp \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
            fprintf(fileID,'recorder Node %s %s/nodal_accel_%s.%s -time -node %i -dof 1 3 accel \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
            fprintf(fileID,'recorder Node %s %s/nodal_reaction_%s.%s -time -node %i -dof 1 3 4 6 reaction \n', file_type, write_dir, num2str(node_ids(i)), file_ext, node_ids(i));
    end
end

% Define Element Recorders
for i = 1:length(element_ids)
    if strcmp(dimension,'2D')
        fprintf(fileID,'recorder Element %s %s/element_force_%s.%s -time -ele %i localForce \n', file_type, write_dir, num2str(element_ids(i)), file_ext, element_ids(i));
    else
        fprintf(fileID,'recorder Element %s %s/element_force_%s.%s -time -ele %i localForce \n', file_type, write_dir, num2str(element_ids(i)), file_ext, element_ids(i));
    end
end


% % Movie Recorders
% fprintf(fileID,'recorder display "Displaced shape" 10 10 500 500 -wipe \n');
% fprintf(fileID,'prp 200.0 50.0 50.0; \n');
% fprintf(fileID,'vup 0.0 1.0 0.0; \n');
% if strcmp(dimension,'2D')
%     fprintf(fileID,'vpn 0.0 0.0 1.0; \n');
% else
%     fprintf(fileID,'vpn 0.4 0.25 1; \n');
% end
% %     fprintf(fileID,'viewWindow -1000 1000 -1000 1000 \n');
% movie_scale = 10;
% fprintf(fileID,'display 1 5 %f \n',movie_scale);


fprintf(fileID,'puts "Defining Recorders Complete" \n');

% Close File
fclose(fileID);

end

