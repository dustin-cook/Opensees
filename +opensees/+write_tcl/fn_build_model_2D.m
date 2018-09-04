function [ node ] = fn_build_model_2D( output_dir, node, element, joint, hinge, analysis, dimension )
%UNTITLED6 Summary of this function goes here

%% Load element properties table
ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

%% Write TCL file
file_name = [output_dir filesep 'model.tcl'];
fileID = fopen(file_name,'w');

%% Define the model (2 dimensions, 3 dof)
if strcmp(dimension,'2D')
    fprintf(fileID,'model basic -ndm 2 -ndf 3 \n');
elseif strcmp(dimension,'3D')
    fprintf(fileID,'model basic -ndm 3 -ndf 6 \n');
else
    error('Number of dimensions not valid')
end

%% Define nodes (inches)
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'node %d %f %f \n',node.id(i),node.x(i),node.y(i));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'node %d %f %f %f \n',node.id(i),node.x(i),node.y(i),node.z(i));
    end
end

%% Set boundary conditions at each node (6dof) (fix = 1, free = 0)
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'fix %d %s %s %s \n',node.id(i),node.fix{i}(2),node.fix{i}(3),node.fix{i}(7));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'fix %d %s %s %s %s %s %s \n',node.id(i),node.fix{i}(2),node.fix{i}(3),node.fix{i}(4),node.fix{i}(5),node.fix{i}(6),node.fix{i}(7));
    end
end

%% Define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'mass %d %f 0. 0. \n',node.id(i), node.mass(i));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'mass %d %f 0. %f 0. 0. 0. \n',node.id(i), node.mass(i), node.mass(i));
    end
end

%% Linear Transformation
if strcmp(dimension,'2D')
    fprintf(fileID,'geomTransf PDelta 1 \n'); % Columns
    fprintf(fileID,'geomTransf PDelta 2 \n'); % Beams (x-direction)
elseif strcmp(dimension,'3D')
    fprintf(fileID,'geomTransf PDelta 1 0 0 1 \n'); % Columns
    fprintf(fileID,'geomTransf PDelta 2 0 0 1 \n'); % Beams (x-direction)
    fprintf(fileID,'geomTransf PDelta 3 -1 0 0 \n'); % Girders (z-direction)
    fprintf(fileID,'geomTransf PDelta 4 1 0 0 \n'); % Columns (z-direction)
end

% Define Elements
for i = 1:height(element)
    % Fundamental element properties
    ele_props = ele_props_table(ele_props_table.id == element.ele_id(i),:);
    
    % Define Geotransform property for this element
    if strcmp(dimension,'2D')
        if strcmp(element.type{i},'column') || strcmp(ele_props.type,'wall')
            geotransf = 1;
        elseif strcmp(element.type{i},'beam')
            geotransf = 2;
        end
    elseif strcmp(dimension,'3D')
        if strcmp(ele_props.type,'column') || strcmp(ele_props.type,'wall')
            if strcmp(element.direction{i},'x')
                geotransf = 1;
            elseif strcmp(element.direction{i},'z')
                geotransf = 4;
            end
        elseif strcmp(ele_props.type,'beam')
            if strcmp(element.direction{i},'x')
                geotransf = 2;
            elseif strcmp(element.direction{i},'z')
                geotransf = 3;
            end
        end
    end
    
    % Beams and Columns Assignment
    if strcmp(element.type{i},'beam') || strcmp(element.type{i},'column') 
        if analysis.nonlinear ~= 0 % Nonlinear Analysis
            Iz_ele = ele_props.iz*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
            % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
            fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,Iz_ele,geotransf);
            
            Iz_ele = ele_props.iz*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
            Iy_ele = ele_props.iy*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod);
            % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
            fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,Iy_ele,Iz_ele,geotransf);
        else
            if analysis.model_type == 1 % SDOF
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),element.a,element.e,element.i,geotransf);
            elseif analysis.model_type == 2 %MDOF
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
            end
        end
    % Wall Assignment
    elseif strcmp(element.type{i},'wall')
        if analysis.nonlinear ~= 0 % EPP Nonlinear Analysis
            %uniaxialMaterial ElasticPP $matTag $E $epsyP
            epsy = ele_props.fc_e/ele_props.e;
            fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(i), ele_props.e, epsy); % Elastic Perfectly Plastic Material
        else
            % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(i),ele_props.e*0.5);
        end
        % section Fiber $secTag <-GJ $GJ> {
        fprintf(fileID,'section Fiber %i { \n',element.id(i));
            % patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
            fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n',element.id(i),10,10,-ele_props.d/2,-ele_props.w/2,ele_props.d/2,ele_props.w/2);
        fprintf(fileID,'} \n');
        
        % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>
        fprintf(fileID,'element forceBeamColumn %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),5,element.id(i),1);  

%         % nDMaterial ElasticIsotropic $matTag $E $v <$rho>
%         fprintf(fileID,'nDMaterial ElasticIsotropic %i %f %f \n',element.id(i),ele_props.e,ele_props.poisson_ratio);
%         % section PlateFiber $secTag $matTag $h
%         fprintf(fileID,'section PlateFiber %i %i %f \n',element.id(i),element.id(i),12);
%         % element ShellMITC4 $eleTag $iNode $jNode $kNode $lNode $secTag
%         fprintf(fileID,'element ShellMITC4 %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),element.node_3(i),element.node_4(i),element.id(i));
    
    % Truss Assigment
    elseif strcmp(element.type{i},'truss')
        % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
        fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(i),ele_props.e);
        % element truss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        fprintf(fileID,'element truss %i %i %i %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,element.id(i));
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
if exist('joint','var')
    for i = 1:height(joint)
        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+1,joint.x_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 2 \n',joint.id(i)*10+2,joint.center(i),joint.x_pos(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+3,joint.y_neg(i),joint.center(i));
        fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 200000. 1 \n',joint.id(i)*10+4,joint.center(i),joint.y_pos(i));
    end
end

% Define Plastic Hinges
if exist('hinge','var')
    for i = 1:height(hinge)
        element.id(end + 1) = element.id(end) + 1;
        ele = element(element.id == hinge.element_id(i),:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        % Hinge Stiffnes Calc from appendix B of Ibarra and Krawinkler 2005
        k_mem = 6*(ele_props.e*ele_props.iz)/ele.length; 
        k_ele = ((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod)*k_mem; 
        k_spring = analysis.hinge_stiff_mod*k_ele;
        if analysis.nonlinear == 1 % IMK Rotational Hinge
            % Load linear element table
            ele_lin_table = readtable([output_dir filesep 'element_linear.csv'],'ReadVariableNames',true);
            ele_lin = ele_lin_table(ele_lin_table.id == hinge.element_id(i),:);
            theta_pc = (ele_lin.b_hinge - ele_lin.a_hinge)/2;
            theta_u_pos = ele_lin.Mn_aci_pos/k_spring + ele_lin.b_hinge;
            theta_u_neg = ele_lin.Mn_aci_neg/k_spring + ele_lin.b_hinge;
            if ele_lin.a_hinge > 0
                as_mem_pos = min(((ele_lin.Mp_pos-ele_lin.Mn_aci_pos)/ele_lin.a_hinge)/k_mem,0.1); % No more than 10% of the elastic stiffness acording to ASCE 41-17 10.3.1.2
                as_sping_pos = ((n+1)*as_mem_pos)/(n+1-n*as_mem_pos);
                as_mem_neg = min(((ele_lin.Mp_neg-ele_lin.Mn_aci_neg)/ele_lin.a_hinge)/k_mem,0.1); % No more than 10% of the elastic stiffness acording to ASCE 41-17 10.3.1.2
                as_sping_neg = ((n+1)*as_mem_neg)/(n+1-n*as_mem_neg);
            else
                as_sping_pos = 0.0;
                as_sping_neg = 0.0;
            end
            %uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
            fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',element.id(end), k_spring, as_sping_pos, as_sping_neg, ele_lin.Mn_aci_pos, -ele_lin.Mn_aci_neg, ele_lin.a_hinge, ele_lin.a_hinge, theta_pc, theta_pc, ele_lin.c_hinge, ele_lin.c_hinge, theta_u_pos, theta_u_neg);
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)); % Element Id for Hinge
            fprintf(fileID,'equalDOF %i %i 1 2 \n',hinge.node_1(i),hinge.node_2(i));
        elseif analysis.nonlinear == 2 % Elastic Perfectly Plastic Rotational Hinges
            epsy = 6.5*min([ele.Mn_aci_pos,ele.Mn_aci_neg])/k_spring;
            %uniaxialMaterial ElasticPP $matTag $E $epsyP
            fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(end), k_spring, epsy); % Elastic Perfectly Plastic Material
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)); % Element Id for Hinge
            fprintf(fileID,'equalDOF %i %i 1 2 \n',hinge.node_1(i),hinge.node_2(i));
        end
    end
end

% Print model to file 
% fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

