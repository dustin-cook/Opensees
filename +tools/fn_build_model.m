function [ node ] = fn_build_model( output_dir, node, element, story, joint )
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

% Linear Transformation
fprintf(fileID,'geomTransf PDelta 1 0 0 1 \n'); % Columns
fprintf(fileID,'geomTransf PDelta 2 0 0 1 \n'); % Beams (x-direction)
fprintf(fileID,'geomTransf PDelta 3 -1 0 0 \n'); % Girders (y-direction)
fprintf(fileID,'geomTransf PDelta 4 1 0 0 \n'); % Columns (y-direction)

% Define Elements (columns and beam)
% element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
for i = 1:length(element.id)
    fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_start(i),element.node_end(i),element.a(i),element.e(i),element.g(i),element.j(i),element.iy(i),element.iz(i),element.orientation(i));
end

% % Define Materials
% %uniaxialMaterial Elastic $matTag $E
% fprintf(fileID,'uniaxialMaterial Elastic 1 999999999. \n'); %Rigid Elastic Material
% 
% % Define Joints
% % element Joint3D %tag %Nx- %Nx+ %Ny- %Ny+ %Nz- %Nz+ %Nc %MatX %MatY %MatZ %LrgDspTag
% for i = 1:length(joint.id)
%     fprintf(fileID,'element Joint3D %i %i %i %i %i %i %i %i 1 1 1 0 \n',joint.id(i),joint.x_neg(i),joint.x_pos(i),joint.y_neg(i),joint.y_pos(i),joint.z_neg(i),joint.z_pos(i),joint.center(i));
% end

% Define Joints as rigid beam-column elements
for i = 1:length(joint.id)
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 2 \n',joint.id(i)*10+1,joint.x_neg(i),joint.center(i));
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 2 \n',joint.id(i)*10+2,joint.center(i),joint.x_pos(i));
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 1 \n',joint.id(i)*10+3,joint.y_neg(i),joint.center(i));
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 1 \n',joint.id(i)*10+4,joint.center(i),joint.y_pos(i));
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 3 \n',joint.id(i)*10+5,joint.z_neg(i),joint.center(i));
    fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999. 99999. 99999. 200000. 200000. 3 \n',joint.id(i)*10+6,joint.center(i),joint.z_pos(i));
end

% Set Node weight to zero
node.weight = zeros(1,length(node.id));

% Define Rigid Slabs and assign weight to slab
for i = 1:length(story.id)
    slab_ht = story.y_offset(i) + story.story_ht(i);
    node_filter = (node.y == slab_ht);% & (node.z ~= 450) & (node.x == 0);
    nodes_on_slab{i} = node.id(node_filter);
    
    % Find slab extreme values
    max_x = max(node.x(node_filter));
    min_x = min(node.x(node_filter));
    max_z = max(node.z(node_filter));
    min_z = min(node.z(node_filter));
    
    % assign each node trib area ratio
    for j = 1:length(nodes_on_slab{i})
        n_id = nodes_on_slab{i}(j);
        if (node.x(n_id) == max_x || node.x(n_id) == min_x) && (node.z(n_id) == max_z || node.z(n_id) == min_z) % Corner Slab Node
            node.trib_area_ratio(n_id) = 0.25;
        elseif (node.x(n_id) == max_x || node.x(n_id) == min_x) || (node.z(n_id) == max_z || node.z(n_id) == min_z) % Side Slab Node
            node.trib_area_ratio(n_id) = 0.5;
        else % Center Slab Node
            node.trib_area_ratio(n_id) = 1;
        end
    end
    
    % total trib ratio
    wt_per_node = story.story_wt(i)/sum(node.trib_area_ratio(node_filter));
    node.weight(node_filter) = node.trib_area_ratio(node_filter) * wt_per_node;

end

node.mass = node.weight*10/386;

% define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
    fprintf(fileID,'mass %d %f 0. %f 0. 0. 0. \n',node.id(i), node.mass(i), node.mass(i));
end

for i = 1:length(story.id)
    % Define Rigid Slabs
    fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(nodes_on_slab{i}));
end
    
% Print model to file 
fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

