function [ node ] = fn_build_model_3D( output_dir, node, element, story, joint, hinge, analysis )
%UNTITLED6 Summary of this function goes here

%% Load element properties table
ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

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
    fprintf(fileID,'fix %d %d %d %d %d %d %d \n',node.id(i),node.fix{i}(1),node.fix{i}(2),node.fix{i}(3),node.fix{i}(4),node.fix{i}(5),node.fix{i}(6));
end

% define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
    fprintf(fileID,'mass %d %f 0. %f 0. 0. 0. \n',node.id(i), node.mass(i), node.mass(i));
end

% Linear Transformation
fprintf(fileID,'geomTransf PDelta 1 0 0 1 \n'); % Columns
fprintf(fileID,'geomTransf PDelta 2 0 0 1 \n'); % Beams (x-direction)
fprintf(fileID,'geomTransf PDelta 3 -1 0 0 \n'); % Girders (y-direction)
fprintf(fileID,'geomTransf PDelta 4 1 0 0 \n'); % Columns (y-direction)

% Define Elements
for i = 1:length(element.id)
    ele_props = ele_props_table(ele_props_table.id == element.ele_id(i),:);
    % Beams, Columns and Rigid Links
    if strcmp(ele_props.type,'beam') || strcmp(ele_props.type,'column') || strcmp(ele_props.type,'rigid link') 
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,element.orientation(i));
    % Walls
    elseif strcmp(ele_props.type,'wall')
        %section ElasticMembranePlateSection $secTag $E $nu $h $rho
        fprintf(fileID,'section ElasticMembranePlateSection %i %f %f %f 0.0 \n',i,ele_props.e,ele_props.poisson_ratio,ele_props.w); %Elastic Wall Section
        %element ShellMITC4 $eleTag $iNode $jNode $kNode $lNode $secTag
        fprintf(fileID,'element ShellMITC4 %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),element.node_3(i),element.node_4(i),i); % Model Wall as shell
    end
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
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*10+1,joint.x_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*10+2,joint.center(i),joint.x_pos(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*10+3,joint.y_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*10+4,joint.center(i),joint.y_pos(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*10+5,joint.z_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*10+6,joint.center(i),joint.z_pos(i));
    end
end

% Define Rigid Slabs
for i = 1:length(story.id)
    fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(story.nodes_on_slab{i}));
end

% Define Plastic Hinges
if isfield(hinge,'id')
    %uniaxialMaterial ElasticPP $matTag $E $epsyP
    fprintf(fileID,'uniaxialMaterial ElasticPP 1 100. 0.5 \n'); % Elastic Perfectly Plastic Material
    %uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
    fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented 2 100. 1. 1. 0. 0. 100. 100. 100. 100. 1. 1. 1. 1. 0.05 0.05 0.05 0.05 0. 0. 0.05 0.05 1. 1. \n');
        
    for i = 1:length(hinge.id)
        element.id(end + 1) = element.id(end) + 1;
        %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
        if analysis.nonlinear == 1 % Shear Spring
            fprintf(fileID,'element zeroLength %i %i %i -mat 1 -dir 1 \n',element.id(end),hinge.node_1(i),hinge.node_2(i)); % Element Id for Hinge
        elseif analysis.nonlinear == 2 % Rotational Spring
            fprintf(fileID,'element zeroLength %i %i %i -mat 2 -dir 6 \n',element.id(end),hinge.node_1(i),hinge.node_2(i)); % Element Id for Hinge
        end
    end
end

% Print model to file 
fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

