function [ node ] = fn_build_model_2D( output_dir, node, element, story, joint, hinge, analysis )
%UNTITLED6 Summary of this function goes here

%% Load element properties table
ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

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
for i = 1:length(element.id)
    ele_props = ele_props_table(ele_props_table.id == element.ele_id(i),:);
    % Beams, Columns and Rigid Links
    if strcmp(ele_props.type,'beam') || strcmp(ele_props.type,'column') || strcmp(ele_props.type,'rigid link') 
        if analysis.nonlinear ~= 0 % Nonlinear Analysis
            ele_props.iz = ele_props.iz*1.2; % Add stiffness to element to account for two springs
        end
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,element.orientation(i));
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
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+1,joint.x_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+2,joint.center(i),joint.x_pos(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+3,joint.y_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+4,joint.center(i),joint.y_pos(i));
    end
end

% Define Plastic Hinges
if isfield(hinge,'id')
% Load linear element table
ele_lin_table = readtable([output_dir filesep 'element_linear.csv'],'ReadVariableNames',true);
    if analysis.nonlinear == 1 % Rotational Spring
        for i = 1:length(hinge.id)
            element.id(end + 1) = element.id(end) + 1;
            ele = element(element.id == hinge.element_id(i),:);
            ele_lin = ele_lin_table(ele_lin_table.id == hinge.element_id(i),:);
            ele_props = ele_props_table(ele_props_table.id == ele_lin.ele_id,:);
            k = 48*(ele_props.e*ele_props.iz)/ele.length;
            theta_pc = (ele_lin.b_hinge - ele_lin.a_hinge)/2;
            theta_u = ele_lin.Mn_aci_c/k + ele_lin.b_hinge;
            strain_hard_ratio = 0.0; %ele_lin.Mp_c/ele_lin.Mn_aci_c;
            %uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
            fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 1000. 1000. 1000. 1000. 1. 1. 1. 1. %f %f %f %f %f %f %f %f 1.0 1.0 \n',element.id(end), k, strain_hard_ratio, strain_hard_ratio, ele_lin.Mn_aci_c, ele_lin.Mn_aci_c, ele_lin.a_hinge, ele_lin.a_hinge, theta_pc, theta_pc, ele_lin.c_hinge, ele_lin.c_hinge, theta_u, theta_u);
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)); % Element Id for Hinge
            fprintf(fileID,'equalDOF %i %i 1 2 \n',hinge.node_1(i),hinge.node_2(i));
        end
    elseif analysis.nonlinear == 2 % Elastic Perfectly Plastic Rotational Hinges
        for i = 1:length(hinge.id)
            element.id(end + 1) = element.id(end) + 1;
            %uniaxialMaterial ElasticPP $matTag $E $epsyP
            fprintf(fileID,'uniaxialMaterial ElasticPP %i 29000000000. 0.001 \n', element.id(end)); % Elastic Perfectly Plastic Material
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)); % Element Id for Hinge
            fprintf(fileID,'equalDOF %i %i 1 2 \n',hinge.node_1(i),hinge.node_2(i));
        end
    end
end

% Print model to file 
fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

