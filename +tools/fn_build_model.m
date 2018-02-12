function [ ] = fn_build_model( output_dir, node, element, story )
%UNTITLED6 Summary of this function goes here

%% Write TCL file
file_name = [output_dir filesep 'model.tcl'];
fileID = fopen(file_name,'w');

% Define the model (3 dimensions, 6 dof)
fprintf(fileID,'model basic -ndm 3 -ndf 6 \n');

% define nodes (inches)
for i = 1:length(node.id)
    fprintf(fileID,'node %d %f %f %f \n',node.id(i),node.x(i),node.y(i),node.z(i));
end

% set boundary conditions at each node (6dof) (fix = 1, free = 0)
for i = 1:length(node.id)
    fprintf(fileID,'fix %d %d %d %d %d %d %d \n',node.id(i),node.fix(i,1),node.fix(i,1),node.fix(i,1),node.fix(i,1),node.fix(i,1),node.fix(i,1));
end

% define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
fprintf(fileID,'mass %d %f 0. %f 0. 0. 0. \n',node.id(i), node.mass(i), node.mass(i));
end

% Linear Transformation
fprintf(fileID,'geomTransf PDelta 1 0 0 1 \n'); % Columns
fprintf(fileID,'geomTransf PDelta 2 0 0 1 \n'); % Beams (x-direction)
fprintf(fileID,'geomTransf PDelta 3 -1 0 0 \n'); % Girders (y-direction)

% Define Elements (columns and beam)
% element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
for i = 1:length(element.id)
    fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_start(i),element.node_end(i),element.a(i),element.e(i),element.g(i),element.j(i),element.iy(i),element.iz(i),element.orientation(i));
end

% Define Rigid Slabs
for i = 1:length(story.id)
    slab_ht = story.y_offset(i) + story.story_ht(i);
    node_filter = (node.y == slab_ht);% & (node.z ~= 450) & (node.x == 0);
    nodes_on_slab = node.id(node_filter);
    fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(nodes_on_slab));
%     fprintf(fileID,'node %d 750 %d 450 \n',9990+i,slab_ht);
%     fprintf(fileID,'rigidDiaphragm 2 %d %s \n',9990+i,num2str([13 41 69 97]));
end

% Print model to file 
fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

