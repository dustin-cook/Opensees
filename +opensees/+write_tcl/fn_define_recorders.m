function [ ] = fn_define_recorders( write_dir, dimension, nodes, element, joint, hinge, analysis )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Define File Type to Write to
if analysis.write_xml
    file_type = '-xml';
    file_ext = 'xml';
else
    file_type = '-file';
    file_ext = 'txt';
end

% Write Recorder File
file_name = [write_dir, filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

%% Dynamic Recorders
if analysis.type == 1 
    % Define Node recorders
    fprintf(fileID,'recorder Node %s %s/nodal_disp_x.%s -time -node %s -dof 1 disp \n', file_type, write_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_accel_x.%s -time -node %s -dof 1 accel \n', file_type, write_dir, file_ext, num2str(nodes));
%     fprintf(fileID,'recorder Node %s %s/nodal_disp_y.%s -time -node %s -dof 2 disp \n', file_type, write_dir, file_ext, num2str(nodes));
%     fprintf(fileID,'recorder Node %s %s/nodal_accel_y.%s -time -node %s -dof 2 accel \n', file_type, write_dir, file_ext, num2str(nodes));
    if strcmp(dimension,'3D')
        fprintf(fileID,'recorder Node %s %s/nodal_disp_z.%s -time -node %s -dof 3 disp \n', file_type, write_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_accel_z.%s -time -node %s -dof 3 accel \n', file_type, write_dir, file_ext, num2str(nodes));
    end

    % Define Element Recorders
    % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
    if strcmp(dimension,'2D')
        fprintf(fileID,'recorder Element %s %s/element_force.%s -time -ele %s -dof 1 2 3 6 localForce \n', file_type, write_dir, file_ext, num2str(element.id'));
    else
        fprintf(fileID,'recorder Element %s %s/element_force.%s -time -ele %s -dof 1 2 6 12 localForce \n', file_type, write_dir, file_ext, num2str(element.id'));
    end
    
    % Hinges
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        fprintf(fileID,'recorder Element %s %s/hinge_force_all.%s -time -ele %s -dof 1 3 4 6 force \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
        fprintf(fileID,'recorder Element %s %s/hinge_deformation_all.%s -time -ele %s deformation \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
    end
    
%% Pushover Recorders
elseif analysis.type == 2 || analysis.type == 3 
    % Nodal Reaction Recorders
    if strcmp(analysis.pushover_direction,'x')
        fprintf(fileID,'recorder Node %s %s/nodal_disp_x.%s -time -node %s -dof 1 disp \n', file_type, write_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_reaction_x.%s -time -node %s -dof 1 reaction \n', file_type, write_dir, file_ext, num2str(nodes));
        if strcmp(dimension,'2D')
            fprintf(fileID,'recorder Element %s %s/element_force_x.%s -time -ele %s -dof 1 2 3 6 localForce \n', file_type, write_dir, file_ext, num2str(element.id'));
        else
            fprintf(fileID,'recorder Element %s %s/element_force_x.%s -time -ele %s -dof 1 2 6 12 localForce \n', file_type, write_dir, file_ext, num2str(element.id'));
        end
    elseif strcmp(analysis.pushover_direction,'z')
        fprintf(fileID,'recorder Node %s %s/nodal_disp_z.%s -time -node %s -dof 3 disp \n', file_type, write_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_reaction_z.%s -time -node %s -dof 3 reaction \n', file_type, write_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Element %s %s/element_force_z.%s -time -ele %s -dof 1 2 6 12 localForce \n', file_type, write_dir, file_ext, num2str(element.id'));
    end
    
    % Hinges
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        if strcmp(analysis.pushover_direction,'x')
            fprintf(fileID,'recorder Element %s %s/hinge_force_x.%s -time -ele %s -dof 1 6 force \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
            fprintf(fileID,'recorder Element %s %s/hinge_deformation_x.%s -time -ele %s deformation \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
        elseif strcmp(analysis.pushover_direction,'z')
            fprintf(fileID,'recorder Element %s %s/hinge_force_z.%s -time -ele %s -dof 3 4 force \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
            fprintf(fileID,'recorder Element %s %s/hinge_deformation_z.%s -time -ele %s deformation \n', file_type, write_dir, file_ext, num2str(element.id(end) + hinge.id'));
        end
    end
end

%% Joints
% if ~isempty(joint)
%     fprintf(fileID,'recorder Element %s %s/joint_force_all.%s -time -ele %s force \n', file_type, write_dir, file_ext, num2str(joint.id' + 10000));
%     fprintf(fileID,'recorder Element %s %s/joint_deformation_all.%s -time -ele %s deformation \n', file_type, write_dir, file_ext, num2str(joint.id' + 10000));
% end

%% Movie Recorders
if analysis.play_movie
    fprintf(fileID,'recorder display "Displaced shape" 10 10 500 500 -wipe \n');
    fprintf(fileID,'prp 200.0 50.0 50.0; \n');
    fprintf(fileID,'vup 0.0 1.0 0.0; \n');
    if strcmp(dimension,'2D')
        fprintf(fileID,'vpn 0.0 0.0 1.0; \n');
    else
        fprintf(fileID,'vpn 0.4 0.25 1; \n');
    end
    %     fprintf(fileID,'viewWindow -1000 1000 -1000 1000 \n');
    fprintf(fileID,'display 1 5 %f \n',analysis.movie_scale);
end

% Close File
fclose(fileID);

end

