function [ node ] = fn_build_model_2D( output_dir, node, element, story, joint, wall, hinge, analysis )
%UNTITLED6 Summary of this function goes here

%% Write TCL file
file_name = [output_dir filesep 'model.tcl'];
fileID = fopen(file_name,'w');

% Define the model (2 dimensions, 3 dof)
fprintf(fileID,'model basic -ndm 2 -ndf 3 \n');

% define nodes (inches)
for i = 1:length(node.id)
    fprintf(fileID,'node %d %f %f \n',node.id(i),node.x(i),node.y(i));
end

% set boundary conditions at each node (6dof) (fix = 1, free = 0)
for i = 1:length(node.id)
    fprintf(fileID,'fix %d %d %d %d \n',node.id(i),node.fix{i}(1),node.fix{i}(2),node.fix{i}(6));
end

% define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
    fprintf(fileID,'mass %d %f 0. 0. \n',node.id(i), node.mass(i));
end

% Linear Transformation
fprintf(fileID,'geomTransf PDelta 1 \n'); % Columns
fprintf(fileID,'geomTransf PDelta 2 \n'); % Beams (x-direction)

% Define Elements (columns and beam)
% element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
if sum(strcmp('id',element.Properties.VariableNames)) > 0
    for i = 1:length(element.id)
        fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_start(i),element.node_end(i),element.a(i),element.e(i),element.iz(i),element.orientation(i));
    end
else
    element.id = 0;
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
if isfield(joint,'id')
    for i = 1:length(joint.id)
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+1,joint.x_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+2,joint.center(i),joint.x_pos(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+3,joint.y_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+4,joint.center(i),joint.y_pos(i));
    end
end

% % Define Rigid Slabs
% for i = 1:length(story.id)
%     fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(story.nodes_on_slab{i}));
% end

% Define Walls Sections
mat_id = 0;
if isfield(wall,'id')
    for i = 1:length(wall.id)
        element.id(end + 1) = element.id(end) + 1;
        mat_id = mat_id + 1;
        %uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg
        fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',mat_id,wall.e(i)); %Elastic Wall Section
        %element quad $eleTag $iNode $jNode $kNode $lNode $thick $type $matTag <$pressure $rho $b1 $b2>
        fprintf(fileID,'element quad %i %i %i %i %i %f PlaneStress %i \n',element.id(end),wall.node_1(i),wall.node_2(i),wall.node_3(i),wall.node_4(i),wall.thickness(i),mat_id); % Model Wall as shell
    end
end

% Define Plastic Hinges
if isfield(hinge,'id')
    mat_id = mat_id + 1;
    if analysis.nonlinear == 1 % Shear Spring
        %uniaxialMaterial ElasticPP $matTag $E $epsyP
        fprintf(fileID,'uniaxialMaterial ElasticPP %i 100. 0.5 \n', mat_id); % Elastic Perfectly Plastic Material
        for i = 1:length(hinge.id)
            element.id(end + 1) = element.id(end) + 1;
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 1 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), mat_id); % Element Id for Hinge
        end
    elseif analysis.nonlinear == 2 % Rotational Spring
        %uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
        fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i 100. 1. 1. 0. 0. 100. 100. 100. 100. 1. 1. 1. 1. 0.05. 0.05. 0.05. 0.05. 0. 0. 0.05 0.05 1. 1. \n',mat_id);
        for i = 1:length(hinge.id)
            element.id(end + 1) = element.id(end) + 1;
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), mat_id); % Element Id for Hinge
        end
    end
end

% Print model to file 
fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

