function [ ] = fn_define_recorders( output_dir, dimension, nodes, element, hinge, analysis )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Define File Type to Write to
if analysis.summit_SP
    file_type = '-xml';
    file_ext = 'xml';
else
    file_type = '-file';
    file_ext = 'txt';
end

% Write Recorder File
file_name = [output_dir, filesep 'recorders.tcl'];
fileID = fopen(file_name,'w');

if analysis.full_recorders == 1
    %% Define Node recorders
    fprintf(fileID,'recorder Node %s %s/nodal_disp_x.%s -time -node %s -dof 1 disp \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_disp_y.%s -time -node %s -dof 2 disp \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_reaction_x.%s -time -node %s -dof 1 reaction \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_reaction_y.%s -time -node %s -dof 2 reaction \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_accel_x.%s -time -node %s -dof 1 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_accel_y.%s -time -node %s -dof 2 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    if strcmp(dimension,'3D')
        fprintf(fileID,'recorder Node %s %s/nodal_disp_z.%s -time -node %s -dof 3 disp \n', file_type, output_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_reaction_z.%s -time -node %s -dof 3 reaction \n', file_type, output_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_accel_z.%s -time -node %s -dof 3 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    end

    %% Define Element Recorders
    % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
    fprintf(fileID,'recorder Element %s %s/element_force.%s -time -ele %s localForce \n', file_type, output_dir, file_ext, num2str(element.id'));

else % Default simpler recorder set
    %% Define Node recorders
    fprintf(fileID,'recorder Node %s %s/nodal_disp_x.%s -time -node %s -dof 1 disp \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_disp_y.%s -time -node %s -dof 2 disp \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_accel_x.%s -time -node %s -dof 1 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    fprintf(fileID,'recorder Node %s %s/nodal_accel_y.%s -time -node %s -dof 2 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    if strcmp(dimension,'3D')
        fprintf(fileID,'recorder Node %s %s/nodal_disp_z.%s -time -node %s -dof 3 disp \n', file_type, output_dir, file_ext, num2str(nodes));
        fprintf(fileID,'recorder Node %s %s/nodal_accel_z.%s -time -node %s -dof 3 accel \n', file_type, output_dir, file_ext, num2str(nodes));
    end

    %% Define Element Recorders
    % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
    if strcmp(dimension,'2D')
        fprintf(fileID,'recorder Element %s %s/element_force.%s -time -ele %s -dof 1 2 3 6 localForce \n', file_type, output_dir, file_ext, num2str(element.id'));
    else
        fprintf(fileID,'recorder Element %s %s/element_force.%s -time -ele %s -dof 1 2 6 12 localForce \n', file_type, output_dir, file_ext, num2str(element.id'));
    end
    
    if analysis.type == 2 || analysis.type == 3 % Default Pushover Recorders
    %% Nodal Reaction Recorders
        if strcmp(analysis.pushover_direction,'x')
    %         fprintf(fileID,'recorder Node -file %s/nodal_disp_x.txt -time -node %s -dof 1 disp \n', output_dir, num2str(nodes));
            fprintf(fileID,'recorder Node %s %s/nodal_reaction_x.%s -time -node %s -dof 1 reaction \n', file_type, output_dir, file_ext, num2str(nodes));
        elseif strcmp(analysis.pushover_direction,'z')
%             fprintf(fileID,'recorder Node -file %s/nodal_disp_z.txt -time -node %s -dof 3 disp \n', output_dir, num2str(nodes));
            fprintf(fileID,'recorder Node %s %s/nodal_reaction_z.%s -time -node %s -dof 3 reaction \n', file_type, output_dir, file_ext, num2str(nodes));
        end
    end
end
% % Walls
% wall_elements = element.id(strcmp(element.type,'wall'))';
% if ~isempty(wall_elements)
%     % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
% %     fprintf(fileID,'recorder Element -file %s/element_force_%d.txt -ele %d Force \n', output_dir, 'walls', num2str(wall_elements));
%     fprintf(fileID,'recorder Element -file %s/element_force_%s.txt -ele %s localForce \n', output_dir, 'walls', num2str(wall_elements));
% end
% 
% % Beams and Columns
% not_wall_elements =  element.id(~strcmp(element.type,'wall'))';
% if ~isempty(not_wall_elements)
%     % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
%     fprintf(fileID,'recorder Element -file %s/element_force_%s.txt -ele %s localForce \n', output_dir, 'beams_and_columns', num2str(not_wall_elements));
% end

% Hinges
if analysis.nonlinear ~= 0 && ~isempty(hinge)
    % recorder Element <-file $fileName> <-time> <-ele ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-region $regTag> <-ele all> ($arg1 $arg2 ...)
    
    % Wall Hinges (shear)
%     wall_elements = element.id(strcmp(element.type,'wall'));
%     shear_hinges = hinge.id(strcmp(hinge.type,'shear'));
%     wall_hinge_ele_ids = element.id(end) + shear_hinges;

%     if ~isempty(wall_hinge_ele_ids)
        fprintf(fileID,'recorder Element %s %s/hinge_force_all.%s -time -ele %s -dof 1 6 force \n', file_type, output_dir, file_ext, num2str(element.id(end) + hinge.id'));
        fprintf(fileID,'recorder Element %s %s/hinge_deformation_all.%s -time -ele %s -dof 1 6 deformation \n', file_type, output_dir, file_ext, num2str(element.id(end) + hinge.id'));
%     end
    
%     % Other Hinges
%     rot_hinges = hinge.id(strcmp(hinge.type,'rotational'));
%     rot_hinge_ele_ids = element.id(end) + rot_hinges;
%     
%     if ~isempty(rot_hinge_ele_ids)
%         fprintf(fileID,'recorder Element -file %s/hinge_moment_all.txt -time -ele %s -dof 6 force \n', output_dir, num2str(rot_hinge_ele_ids));
%         fprintf(fileID,'recorder Element -file %s/hinge_rotation_all.txt -time -ele %s -dof 6 deformation \n', output_dir, num2str(rot_hinge_ele_ids));
%     end
end

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

