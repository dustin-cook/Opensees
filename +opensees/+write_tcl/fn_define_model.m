function [ node ] = fn_define_model( write_dir, node, element, joint, hinge, analysis, dimension, story, read_dir_analysis )
%UNTITLED6 Summary of this function goes here

%% Import Tools
import asce_41.*

%% Load element properties table
ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

%% Write TCL file
file_name = [write_dir filesep 'model.tcl'];
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
        fprintf(fileID,'node %i %f %f \n',node.id(i),node.x(i),node.y(i));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'node %i %f %f %f \n',node.id(i),node.x(i),node.y(i),node.z(i));
    end
end

%% Set boundary conditions at each node (6dof) (fix = 1, free = 0)
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'fix %i %s %s %s \n',node.id(i),node.fix{i}(2),node.fix{i}(3),node.fix{i}(7));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'fix %i %s %s %s %s %s %s \n',node.id(i),node.fix{i}(2),node.fix{i}(3),node.fix{i}(4),node.fix{i}(5),node.fix{i}(6),node.fix{i}(7));
    end
end

%% Define nodal masses (horizontal) (k-s2/in)
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'mass %i %f 0. 0. \n',node.id(i), node.mass(i));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'mass %i %f 0. %f 0. 0. 0. \n',node.id(i), node.mass(i), node.mass(i));
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
    
    %% Beams and Columns Assignment
    if strcmp(element.type{i},'beam') || strcmp(element.type{i},'column') 
        if analysis.nonlinear ~= 0 % Nonlinear Analysis
            Iz_ele = ele_props.iz*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
            Iy_ele = ele_props.iy*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod);
            if strcmp(dimension,'2D')
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,Iz_ele,geotransf);
            elseif strcmp(dimension,'3D')
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,Iy_ele,Iz_ele,geotransf);
            end
        else
            if analysis.model_type == 1 % SDOF
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),element.a,element.e,element.i,geotransf);
            elseif analysis.model_type == 2 %MDOF
                if strcmp(dimension,'2D')
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
                elseif strcmp(dimension,'3D')
%                     fprintf(fileID,'element ElasticTimoshenkoBeam %i %i %i %f %f %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.e,ele_props.g,ele_props.a,ele_props.j,ele_props.iy,ele_props.iz,(5/6)*ele_props.a,(5/6)*ele_props.a,geotransf);
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
                end
            end
        end
        
    %% Wall Assignment
    elseif strcmp(element.type{i},'wall')
        % Elastic
        if analysis.nonlinear == 0 
            if strcmp(dimension,'2D')
                 fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
            elseif strcmp(dimension,'3D')
                fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
            end
        % Shear springs
        elseif analysis.nonlinear == 1
            Iz_ele = ele_props.iz;%*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
            if strcmp(dimension,'2D')
                fprintf(fileID,'element elasticBeamColumn %d %d %d %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
            elseif strcmp(dimension,'3D')
                fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
            end
        % Explicit steel and concrete fibers
        elseif analysis.nonlinear == 2 
            % uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
            fprintf(fileID,'uniaxialMaterial Steel02 %i %f %f 0.05 15. 0.925 0.15 \n', 1000 + element.id(i), ele_props.fy_e, ele_props.Es);
            % uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
            fprintf(fileID,'uniaxialMaterial Concrete04 %i %f -0.002 -0.06 %f %f 0.00015 \n', element.id(i), -ele_props.fc_e, ele_props.e/.35, 7.5*sqrt(ele_props.fc_e));
                
            % section Fiber $secTag <-GJ $GJ> {
            fprintf(fileID,'section Fiber %i { \n',element.id(i));
                % patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
                fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n',element.id(i),round(ele_props.d),round(ele_props.w),0,0,ele_props.d,ele_props.w);
                if analysis.nonlinear == 2
                    As = str2double(strsplit(strrep(strrep(ele_props.As{1},'[',''),']','')));
                    num_bars = str2double(strsplit(strrep(strrep(ele_props.n_b{1},'[',''),']','')));
                    depth_bars = str2double(strsplit(strrep(strrep(ele_props.As_d{1},'[',''),']','')));
                    for row = 1:length(As)
                        % layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
                        height_bars = depth_bars(row);
                        width_start = ele_props.clear_cover + 1;
                        width_end = ele_props.w - width_start;
                        fprintf(fileID,'layer straight %i %i %f %f %f %f %f \n', 1000 + element.id(i), num_bars(row), As(row), height_bars, width_start, height_bars, width_end);
                    end
                end

            fprintf(fileID,'} \n');

            % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>
            fprintf(fileID,'element forceBeamColumn %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),5,element.id(i),geotransf);  
        end

    %% Truss Assigment
    elseif strcmp(element.type{i},'truss')
        % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
        fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(i),ele_props.e);
        % element truss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        fprintf(fileID,'element truss %i %i %i %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,element.id(i));
    end
end

%% Define Joints as rigid
if height(joint) > 0
    % Load in joint properties
    if analysis.nonlinear ~= 0 % Nonlinear analysis
        joint_analysis_temp = load([read_dir_analysis filesep 'joint_analysis.mat']);
        joint_analysis = joint_analysis_temp.joint;
    end
    
    fprintf(fileID,'uniaxialMaterial Elastic 1 999999999999999. \n'); % Rigid Elastic Material
    
    % GO through each joint
    for i = 1:height(joint)
        % Define joint material
%         if analysis.nonlinear ~= 0 % Nonlinear 
%             [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', joint_analysis.Mn(i), joint_analysis.Mn(i), inf, inf, joint_analysis.h(i), joint_analysis.e(i), joint_analysis.iz(i), joint_analysis(i,:), NaN, 0.04 );
%             Ko = moment_vec_pos(1)/rot_vec_pos(1);
%             as_sping_pos = (moment_vec_pos(2)-moment_vec_pos(1))/(rot_vec_pos(2)-rot_vec_pos(1))/Ko;
%             as_sping_neg = (moment_vec_neg(2)-moment_vec_neg(1))/(rot_vec_neg(2)-rot_vec_neg(1))/Ko;
%             theta_pc = rot_vec_neg(3) - rot_vec_neg(2) + joint_analysis.c_hinge(i)*(rot_vec_neg(3) - rot_vec_neg(2))/(1-joint_analysis.c_hinge(i)); % theta pc defined all the way to zero where b defined to residual kink
%             % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%             fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',joint.id(i)+10000, Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc, theta_pc, joint_analysis.c_hinge(i), joint_analysis.c_hinge(i), 999, 999); % Keep residual strength forever
%         end
        if strcmp(dimension,'2D')
            if analysis.joint_model == 1 % Elastic beam column elements
           
                joint_center.x = node.x(node.id == joint.y_pos(i),:);
                joint_center.y = node.y(node.id == joint.x_pos(i),:);
                fprintf(fileID,'node %i %f %f \n',40000+i,joint_center.x,joint_center.y);
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d 100000. 9999999999. 20000000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %d %d %d 100000. 9999999999. 20000000. 2 \n',joint.id(i)*1000000+2,40000+i,joint.x_pos(i));
                fprintf(fileID,'element elasticBeamColumn %d %d %d 100000. 9999999999. 20000000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %d %d %d 100000. 9999999999. 20000000. 1 \n',joint.id(i)*1000000+4,40000+i,joint.y_pos(i));
            elseif analysis.joint_model == 2  % Joint 2D
                % element Joint2D $eleTag $Nd1 $Nd2 $Nd3 $Nd4 $NdC <$Mat1 $Mat2 $Mat3 $Mat4> $MatC $LrgDspTag
                fprintf(fileID,'element Joint2D %i %i %i %i %i %i 1 0 \n', 10000+i, joint.y_pos(i), joint.x_neg(i), joint.y_neg(i), joint.x_pos(i), 10000+i);
            end
        elseif strcmp(dimension,'3D')
            if analysis.joint_model == 1 % Elastic beam column elements
                joint_center.x = node.x(node.id == joint.y_pos(i),:);
                joint_center.y = node.y(node.id == joint.x_pos(i),:);
                joint_center.z = node.z(node.id == joint.x_pos(i),:);
                fprintf(fileID,'node %i %f %f %f \n',40000+i,joint_center.x,joint_center.y,joint_center.z);
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+2,40000+i,joint.x_pos(i));
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,40000+i,joint.y_pos(i));
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*1000000+5,joint.z_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*1000000+6,40000+i,joint.z_pos(i));
            elseif analysis.joint_model == 2  % Joint 3D
                if analysis.nonlinear ~= 0 && analysis.joint_explicit % Nonlinear
                    % element Joint3D %tag %Nx- %Nx+ %Ny- %Ny+ %Nz- %Nz+ %Nc %MatX %MatY %MatZ %LrgDspTag
                    fprintf(fileID,'element Joint3D %i %i %i %i %i %i %i %i %i 1 1 0 \n', 10000+i, joint.x_neg(i), joint.x_pos(i), joint.y_neg(i), joint.y_pos(i), joint.z_neg(i), joint.z_pos(i), 10000+i, 10000+i); 
                else % Linear Rigid
                    % element Joint3D %tag %Nx- %Nx+ %Ny- %Ny+ %Nz- %Nz+ %Nc %MatX %MatY %MatZ %LrgDspTag
                    fprintf(fileID,'element Joint3D %i %i %i %i %i %i %i %i 1 1 1 0 \n', 10000+i, joint.x_neg(i), joint.x_pos(i), joint.y_neg(i), joint.y_pos(i), joint.z_neg(i), joint.z_pos(i), 10000+i);
                end
            end
        end
    end
end

%% Define Rigid Slabs
if analysis.rigid_diaphram
    fprintf(fileID,'uniaxialMaterial Elastic 2 999999999999. \n'); % Rigid Elastic Material
    for i = 1:height(story)
        slab_nodes_at_story = node.id(node.on_slab == story.id(i))';
        frame_nodes_at_story = node.id(node.story == story.id(i))';
        fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(slab_nodes_at_story));
        
        % Create zero length element to connect rigid diaphram to nodes
        for j = 1:length(slab_nodes_at_story)
            new_node = slab_nodes_at_story(j);
            old_node = frame_nodes_at_story(j);
            fprintf(fileID,'element zeroLength %i %i %i -mat 2 2 2 2 2 2 -dir 1 2 3 4 5 6 \n',20000+i*1000+j,new_node,old_node); % Element Id for Hinge
        end
    end
end
                
%% Define Plastic Hinges
if height(hinge) > 0
    % Load linear element table
    if analysis.nonlinear ~= 0
        element_analysis_temp = load([read_dir_analysis filesep 'element_analysis.mat']);
        element_analysis = element_analysis_temp.element;
    end
    
    for i = 1:height(hinge)
        element.id(end + 1) = element.id(end) + 1;
        if strcmp(hinge.type(i),'foundation')
            pile_rot_stiff = 881511412050; % Force-Displacement rotational Stiffness of Bundle of Piles
            pile_lat_stiff = 18138095; % Force-Displacement lateral Stiffness of Bundle of Piles
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(end),pile_rot_stiff); % Rigid Elastic Material
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(end)+8000,pile_lat_stiff); % Rigid Elastic Material
            if strcmp(dimension,'3D') % Input as rotational for now
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i %i %i -dir 1 3 4 6 \n',element.id(end), hinge.node_1(i), hinge.node_2(i), element.id(end)+8000, element.id(end)+8000, element.id(end), element.id(end)); % Element Id for Hinge
                fprintf(fileID,'equalDOF %i %i 2 5 \n',hinge.node_2(i),hinge.node_1(i));
            else
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i -dir 1 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)+8000, element.id(end)); % Element Id for Hinge
                fprintf(fileID,'equalDOF %i %i 2 \n',hinge.node_2(i),hinge.node_1(i));
            end
        else
            ele = element(element.id == hinge.element_id(i),:);
            ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
            if strcmp(ele_props.type,'beam') || strcmp(ele_props.type,'column')  % IMK Rotational Hinge
                % Load in backbone curve
                hinge_props = element_analysis(element_analysis.id == hinge.element_id(i),:);
                [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', hinge_props.Mn_pos, hinge_props.Mn_neg, hinge_props.Mp_pos, hinge_props.Mp_neg, ele.length, ele_props.e, ele_props.iz, hinge_props, analysis.hinge_stiff_mod, 0.1 );

                % Define IMK Parameters
                Ko = moment_vec_pos(1)/rot_vec_pos(1);
                as_sping_pos = (moment_vec_pos(2)-moment_vec_pos(1))/(rot_vec_pos(2)-rot_vec_pos(1))/Ko;
                as_sping_neg = (moment_vec_neg(2)-moment_vec_neg(1))/(rot_vec_neg(2)-rot_vec_neg(1))/Ko;
                theta_pc = rot_vec_neg(3) - rot_vec_neg(2) + hinge_props.c_hinge*(rot_vec_neg(3) - rot_vec_neg(2))/(1-hinge_props.c_hinge); % theta pc defined all the way to zero where b defined to residual kink
                % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
                fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',element.id(end), Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc, theta_pc, hinge_props.c_hinge, hinge_props.c_hinge, 0.9999999, 0.9999999); % Keep residual strength forever

                % uniaxialMaterial Bilin              $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%                 fprintf(fileID,'uniaxialMaterial Bilin %i %f %f %f %f %f 1000.0 1000.0 1000.0 1000.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',element.id(end), Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc, theta_pc, hinge_props.c_hinge, hinge_props.c_hinge, 0.9999999, 0.9999999); % Keep residual strength forever

                % uniaxialMaterial MultiLinear $matTag $u1 $f1 $u2 $f2 $u3 $f3 $u4 $f4 ...
%                 fprintf(fileID,'uniaxialMaterial MultiLinear %i %f %f %f %f %f %f %f %f \n',element.id(end),rot_vec_pos(1),moment_vec_pos(1),rot_vec_pos(2),moment_vec_pos(2),rot_vec_pos(3),moment_vec_pos(3), 9999999, moment_vec_pos(3)); % continue hinge at residual strength
                
                % uniaxialMaterial ElasticPP $matTag $E $epsyP <$epsyN $eps0>
%                 yeild_rot = moment_vec_pos(1)/Ko;
%                 fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(end), Ko, yeild_rot);

                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                if strcmp(dimension,'3D') % Input as rotational for now   
                    % Define Nonlinear Hinge in the weak direction
                    % uniaxialMaterial ElasticPP $matTag $E $epsyP <$epsyN $eps0>
                    out_of_plane_K0 = (analysis.hinge_stiff_mod+1)*6*ele_props.e*ele_props.iy/ele.length;
                    yeild_rot = moment_vec_pos(1)/Ko; % use the same yeild rotation as the primary direction
                    fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(end)+100000, out_of_plane_K0, yeild_rot);

                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',element.id(end), hinge.node_1(i), hinge.node_2(i), element.id(end)); % Element Id for Hinge
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',element.id(end)+100000, hinge.node_1(i), hinge.node_2(i), element.id(end)+100000); % Element Id for Hinge
                    fprintf(fileID,'equalDOF %i %i 1 2 3 5 \n', hinge.node_2(i), hinge.node_1(i));
                else
                end
            elseif strcmp(ele_props.type,'wall')
                % Define Stiffness
                elastic_shear_stiffness = ele_props.g*ele_props.av/ele.length;
                    
                if analysis.nonlinear == 0 % Elastic Lateral Spring for shear deformations
                    % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
                    fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(end)+100,elastic_shear_stiffness);
                    
                    %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                    if strcmp(ele.direction,'x')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 1 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)+100); % Element Id for Hinge
                        fprintf(fileID,'equalDOF %i %i 2 3 4 5 6 \n',hinge.node_2(i),hinge.node_1(i));
                    elseif strcmp(ele.direction,'z')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)+100); % Element Id for Hinge
                        fprintf(fileID,'equalDOF %i %i 1 2 4 5 6 \n',hinge.node_2(i),hinge.node_1(i));
                    end
                    
                elseif analysis.nonlinear == 1 % Multilinear shear spring hinge with min/max limit
                    hinge_props = element_analysis(element_analysis.id == hinge.element_id(i),:);

                    % Define backbone coordinates
                    [ force_vec, disp_vec ] = fn_define_backbone_shear( hinge_props.Vn, ele.length, ele_props.g, ele_props.av, hinge_props );

                    % uniaxialMaterial MultiLinear $matTag $u1 $f1 $u2 $f2 $u3 $f3 $u4 $f4 ...
                    fprintf(fileID,'uniaxialMaterial MultiLinear %i %f %f %f %f %f %f %f %f %f %f \n',element.id(end) + 9000,disp_vec(1),force_vec(1),disp_vec(2),force_vec(2),disp_vec(3),force_vec(3),disp_vec(4),force_vec(4), 999, force_vec(4)); % continue hinge at residual strength

                    % uniaxialMaterial MinMax $matTag $otherTag <-min $minStrain> <-max $maxStrain>
%                     fprintf(fileID,'uniaxialMaterial MinMax %i %i -min %f -max %f \n',element.id(end)+90000,element.id(end) + 9000,-999,999); % Only reduce to zero strength at really high  displacements

                    out_of_plane_K0 = (analysis.hinge_stiff_mod+1)*6*ele_props.e*ele_props.iy/ele.length;
                    yeild_rot = force_vec(1)*ele.length/Ko; % use the same yeild rotation as the primary direction, but converted to flexure
                            
                    %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                    if strcmp(ele.direction,'x')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 1 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)+9000); % Element Id for Hinge
                        if strcmp(dimension,'3D') 
                            fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(end)+100000, out_of_plane_K0, yeild_rot);
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',element.id(end)+100000, hinge.node_1(i), hinge.node_2(i), element.id(end)+100000); % Element Id for Hinge
                            fprintf(fileID,'equalDOF %i %i 2 3 5 6 \n',hinge.node_2(i),hinge.node_1(i));
                        end
                    elseif strcmp(ele.direction,'z')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',element.id(end),hinge.node_1(i),hinge.node_2(i), element.id(end)+9000); % Element Id for Hinge
                        if strcmp(dimension,'3D') 
                            fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n', element.id(end)+100000, out_of_plane_K0, yeild_rot);
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',element.id(end)+100000, hinge.node_1(i), hinge.node_2(i), element.id(end)+100000); % Element Id for Hinge
                            fprintf(fileID,'equalDOF %i %i 1 2 4 5 \n',hinge.node_2(i),hinge.node_1(i));
                        end
                    end
                end
            end
        end
    end
end

% Print model to file 
% fprintf(fileID,'print -file %s/model.txt \n',output_dir);

% Close File
fclose(fileID);

end

