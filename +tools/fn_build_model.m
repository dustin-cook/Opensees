function [ ] = fn_build_model( output_dir, node, element )
%UNTITLED6 Summary of this function goes here
% ASSUMPTIONS:
% 2d linear model single story, 1 bay, uniform members


%% Write TCL file
file_name = [output_dir filesep 'model.tcl'];
fileID = fopen(file_name,'w');

% Define the model (2 dimensions, 3 dof
fprintf(fileID,'#  Define the model (2 dimensions, 3 dof \n');
fprintf(fileID,'model basic -ndm 2 -ndf 3 \n');

% define nodes (inches)
fprintf(fileID,'# define nodes (inches) \n');
for i = 1:length(node.id)
    fprintf(fileID,'node %d %f %f \n',node.id(i),node.x_coor(i),node.y_coor(i));
end

% set boundary conditions at each node (3dof) (fix = 1, free = 0)
fprintf(fileID,'# set boundary conditions at each node (3dof) (fix = 1, free = 0) \n');
for i = 1:length(node.id)
    fprintf(fileID,'fix %d %d %d %d \n',node.id(i),node.fix{i}(1),node.fix{i}(2),node.fix{i}(3));
end

% define nodal masses (horizontal) (units?)
fprintf(fileID,'# define nodal masses (horizontal) (units?) \n');
for i = 1:length(node.id)
fprintf(fileID,'mass %d %f 0. 0. \n',node.id(i), node.mass(i));
end

% Linear Transformation
fprintf(fileID,'# Linear Transformation \n');
fprintf(fileID,'geomTransf Linear 1 \n');

% Define Elements (columns and beam)
% element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag>
fprintf(fileID,'# Define Elements (columns and beam) \n');
fprintf(fileID,'# element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag> \n');
for i = 1:length(element.id)
    fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f 1 \n',element.id(i),element.node_start(i),element.node_end(i),element.A(i),element.E(i),element.I(i));
end

% Close File
fclose(fileID);

end

