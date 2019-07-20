function [ joint_ele_ids ] = fn_define_model( write_dir, node, element, joint, hinge, analysis, dimension, story, read_dir_analysis )
%UNTITLED6 Summary of this function goes here

%% Import Tools
import asce_41.*

%% Load element properties table
ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

%% Write TCL file
file_name = [write_dir filesep 'model.tcl'];
fileID = fopen(file_name,'w');

fprintf(fileID,'puts "Building Model..." \n');

%% Define the model (2 dimensions, 3 dof)
if strcmp(dimension,'2D')
    fprintf(fileID,'model basic -ndm 2 -ndf 3 \n');
    rot_dof_x = 3;
elseif strcmp(dimension,'3D')
    fprintf(fileID,'model basic -ndm 3 -ndf 6 \n');
    rot_dof_x = 6;
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
add_mass = 0.5*(150/386); % one cubic foot of concrete
new_mass = max(node.mass,add_mass);
prcnt_chng = (sum(new_mass) - sum(node.mass)) / sum(node.mass);
if prcnt_chng > 0.01
    error('added too much mass')
end
for i = 1:length(node.id)
    if strcmp(dimension,'2D')
        fprintf(fileID,'mass %i %f %f 0. \n',node.id(i), new_mass(i), new_mass(i));
    elseif strcmp(dimension,'3D')
        fprintf(fileID,'mass %i %f %f %f 0. 0. 0. \n',node.id(i), new_mass(i), new_mass(i), new_mass(i));
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
        if analysis.model_type == 1 % SDOF
            % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
            fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),element.a,element.e,element.iz,geotransf);
        elseif analysis.model_type == 2 %MDOF
            if analysis.nonlinear ~= 0 % Nonlinear Analysis
                Iz_ele = ele_props.iz*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
                Iy_ele = ele_props.iy*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod);
                if strcmp(dimension,'2D')
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,Iz_ele,geotransf);
                elseif strcmp(dimension,'3D')
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,Iy_ele,Iz_ele,geotransf);
                end
            else
                if strcmp(dimension,'2D')
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
                elseif strcmp(dimension,'3D')
%                     fprintf(fileID,'element ElasticTimoshenkoBeam %i %i %i %f %f %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.e,ele_props.g,ele_props.a,ele_props.j,ele_props.iy,ele_props.iz,(5/6)*ele_props.a,(5/6)*ele_props.a,geotransf);
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
                end
            end
        end
        
    %% Wall Assignment
    elseif strcmp(element.type{i},'wall')
        % Elastic
        if analysis.nonlinear == 0 
            if strcmp(dimension,'2D')
                 fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
            elseif strcmp(dimension,'3D')
                fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
            end
        elseif analysis.nonlinear == 1 
            if ~analysis.fiber_walls % Shear springs
                Iz_ele = ele_props.iz;%*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
                elseif strcmp(dimension,'3D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
                end
            else % Explicit steel and concrete fibers 
                % uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
                fprintf(fileID,'uniaxialMaterial Steel02 %i %f %f 0.01 15.0 0.925 0.15 \n', 1000 + element.id(i), ele_props.fy_e, ele_props.Es);
                % uniaxialMaterial Steel01 $matTag $Fy $E0 $b 
%                 fprintf(fileID,'uniaxialMaterial Steel01 %i %f %f 0.1 \n', 1000 + element.id(i), ele_props.fy_e, ele_props.Es*0.35);
                % uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
%                 fprintf(fileID,'uniaxialMaterial Concrete04 %i %f -0.002 -0.06 %f %f 0.00015 \n',2000 + element.id(i), -ele_props.fc_e, ele_props.e*0.35, 7.5*sqrt(ele_props.fc_e));
                % uniaxialMaterial Concrete01 $matTag $fpc $epsc0 $fpcu $epsU
%                 fprintf(fileID,'uniaxialMaterial Concrete01 %i %f -0.002 0.0 -0.005 \n', 2000 + element.id(i), -ele_props.fc_e);
                % uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
%                 fprintf(fileID,'uniaxialMaterial Concrete02 %i %f -0.002 0.0 -0.005 0.75 %f %f \n', 2000 + element.id(i), -ele_props.fc_e, 7.5*sqrt(ele_props.fc_e), ele_props.e*0.1);
                fprintf(fileID,'uniaxialMaterial Concrete04 %i %f -0.002 -0.006 %f %f 0.0001 \n',2000 + element.id(i), -ele_props.fc_e, ele_props.e, 7.5*sqrt(ele_props.fc_n));                
                % section Fiber $secTag <-GJ $GJ> {
                fprintf(fileID,'section Fiber %i { \n',element.id(i));
                    % patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
                    As = str2double(strsplit(strrep(strrep(ele_props.As{1},'[',''),']','')));
                    num_bars = str2double(strsplit(strrep(strrep(ele_props.n_b{1},'[',''),']','')));
                    depth_bars = str2double(strsplit(strrep(strrep(ele_props.As_d{1},'[',''),']','')));
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n',2000 + element.id(i),length(As),3,-ele_props.h/2,-ele_props.w/2,ele_props.h/2,ele_props.w/2);
                    
                    % layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
                    row_ht(1) = (ele_props.clear_cover + 1) - ele_props.w/2;
                    row_ht(2) = ele_props.w/2 - (ele_props.clear_cover + 1);
                    for row = 1:num_bars(1)
                        fprintf(fileID,'layer straight %i %i %f %f %f %f %f \n', 1000 + element.id(i), length(As), mean(As), depth_bars(1)-ele_props.h/2, row_ht(row), depth_bars(end)-ele_props.h/2, row_ht(row));
                    end
%                     for row = 1:length(As)
%                         % layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
%                         height_bars = depth_bars(row);
%                         width_start = ele_props.clear_cover + 1;
%                         width_end = ele_props.w - width_start;
%                         fprintf(fileID,'layer straight %i %i %f %f %f %f %f \n', 1000 + element.id(i), num_bars(row), As(row), height_bars, width_start, height_bars, width_end);
%                     end
                fprintf(fileID,'} \n');

                % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>
                fprintf(fileID,'element forceBeamColumn %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),5,element.id(i),geotransf);  
            end
        end

    %% Truss Assigment
    elseif strcmp(element.type{i},'truss')
        % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
        fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(i),ele_props.e);
        % element truss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        fprintf(fileID,'element truss %i %i %i %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,element.id(i));
    end
end

%% Define Joints
fprintf(fileID,'uniaxialMaterial Elastic 1 999999999. \n'); % Rigid Elastic Material
fprintf(fileID,'uniaxialMaterial Elastic 222222 999999999999. \n'); % Rigid Elastic Material
fprintf(fileID,'uniaxialMaterial Elastic 2 999999999999. \n'); % Rigid Elastic Material
joint_ele_ids = [];
if height(joint) > 0
    % Load in joint properties
    joint_file = [read_dir_analysis filesep 'joint_analysis.mat'];
    if exist(joint_file,'file')
        joint_analysis_temp = load(joint_file);
        joint_analysis = joint_analysis_temp.joint;
    end
    
    % GO through each joint
    for i = 1:height(joint)
        % Define joint material
        if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1 % Nonlinear Joints
                                                                                                     % type, Mn_pos,             Mn_neg,               Mp_pos, Mp_neg, length,                       e,                    iz,    a_hinge,                   b_hinge,                  c_hinge,                    n, strain_harden_ratio
            [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', joint_analysis.Mn(i), joint_analysis.Mn(i), inf, inf, joint_analysis.h(i), joint_analysis.e(i), joint_analysis.iz(i), joint_analysis.a_hinge(i), joint_analysis.b_hinge(i), joint_analysis.c_hinge(i), NaN, 0.04, 'shear' );
            K0 = 1000*moment_vec_pos(1)/rot_vec_pos(1);
            as_sping_pos = (moment_vec_pos(2)-moment_vec_pos(1))/(rot_vec_pos(2)-rot_vec_pos(1))/K0;
            as_sping_neg = (moment_vec_neg(2)-moment_vec_neg(1))/(rot_vec_neg(2)-rot_vec_neg(1))/K0;
            theta_pc_pos = rot_vec_pos(3) - rot_vec_pos(2) + joint_analysis.c_hinge(i)*(rot_vec_pos(3) - rot_vec_pos(2))/(1-joint_analysis.c_hinge(i)); % theta pc defined all the way to zero where b defined to residual kink
            theta_pc_neg = rot_vec_neg(3) - rot_vec_neg(2) + joint_analysis.c_hinge(i)*(rot_vec_neg(3) - rot_vec_neg(2))/(1-joint_analysis.c_hinge(i)); % theta pc defined all the way to zero where b defined to residual kink
            % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
            fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',joint.id(i)+10000, K0, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc_pos, theta_pc_neg, joint_analysis.c_hinge(i), joint_analysis.c_hinge(i), 0.999, 0.999); % Keep residual strength forever
%         else % Linear
%             fprintf(fileID,'uniaxialMaterial Elastic %i 999999999999999. \n',joint.id(i)+10000); % Rigid Elastic Material
        end
        if strcmp(dimension,'2D')
            if analysis.joint_model == 1 % Elastic beam column elements
                joint_ele_ids = [(1:4)+joint.id(i)*1000000,joint_ele_ids];
                joint_center.x = node.x(node.id == joint.y_pos(i),:);
                joint_center.y = node.y(node.id == joint.x_pos(i),:);
                fprintf(fileID,'node %i %f %f \n',40000+i,joint_center.x,joint_center.y);
                fprintf(fileID,'node %i %f %f %f \n',50000+i,joint_center.x,joint_center.y,joint_center.z);
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                fprintf(fileID,'element elasticBeamColumn %i %i %i 100000. 999999999. 20000000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),40000+i);
                fprintf(fileID,'element elasticBeamColumn %i %i %i 100000. 999999999. 20000000. 2 \n',joint.id(i)*1000000+2,40000+i,joint.x_pos(i));
                fprintf(fileID,'element elasticBeamColumn %i %i %i 100000. 999999999. 20000000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),50000+i);
                fprintf(fileID,'element elasticBeamColumn %i %i %i 100000. 999999999. 20000000. 1 \n',joint.id(i)*1000000+4,50000+i,joint.y_pos(i));
                % Scissor Hinge
                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n', joint.id(i)+10000, 40000+i, 50000+i, joint.id(i)+10000); % Currently assumes x direction joint
                fprintf(fileID,'equalDOF %i %i 1 2 3 4 5 \n', 40000+i, 50000+i);
            elseif analysis.joint_model == 2  % Joint 2D
                % element Joint2D $eleTag $Nd1 $Nd2 $Nd3 $Nd4 $NdC <$Mat1 $Mat2 $Mat3 $Mat4> $MatC $LrgDspTag
                fprintf(fileID,'element Joint2D %i %i %i %i %i %i 1 0 \n', 10000+i, joint.y_pos(i), joint.x_neg(i), joint.y_neg(i), joint.x_pos(i), 10000+i);
            end
        elseif strcmp(dimension,'3D')
            if analysis.joint_model == 1 % Elastic beam column elements
                joint_ele_ids = [(1:6)+joint.id(i)*1000000,joint_ele_ids];
                joint_center.x = node.x(node.id == joint.y_pos(i),:);
                joint_center.y = node.y(node.id == joint.x_pos(i),:);
                joint_center.z = node.z(node.id == joint.x_pos(i),:);
                if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1 %  Explicit Joint Model
                    joint_center_node_beams = 40000+i;
                    joint_center_node_columns = 50000+i;
                    fprintf(fileID,'node %i %f %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y,joint_center.z);
                    fprintf(fileID,'node %i %f %f %f \n',joint_center_node_columns,joint_center.x,joint_center.y,joint_center.z);
                else
                    joint_center_node_beams = 40000+i;
                    joint_center_node_columns = 40000+i;
                    fprintf(fileID,'node %i %f %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y,joint_center.z);
                end

                % Joint Fixity Condition
                if exist('joint_analysis','var')
                    joint_rigidity = joint_analysis.implicit_stiff(i);
                else
                    joint_rigidity = 3;
                end
                
                % Grab element info
                bm_left = element(element.id == joint.beam_left(i),:);
                bm_right = element(element.id == joint.beam_right(i),:);
                col_low = element(element.id == joint.column_low(i),:);
                col_high = element(element.id == joint.column_high(i),:);
                
                % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                % Primary Beam Offsets
                if joint_rigidity == 1 % Column Offsets Rigid
                    % Beam Offsets
                    if ~isempty(bm_left)
                        ele_props = ele_props_table(ele_props_table.id == bm_left.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                    end
                    if ~isempty(bm_right)
                        ele_props = ele_props_table(ele_props_table.id == bm_right.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                    end
                    % Column Offsets
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                elseif joint_rigidity == 2 % Beam Offsets Rigid
                    % Primary Beam Offsets
                    if ~isempty(bm_left)
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams);
                    end
                    if ~isempty(bm_right)
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i));
                    end
                    % Column Offsets
                    if ~isempty(col_low)
                        ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    end
                    if ~isempty(col_high)
                        ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                    end
                elseif joint_rigidity == 3 % Both Offsets are Half Rigid
                    % Primary Beam Offsets
                    % Beam Left
                    if ~isempty(bm_left)
                        joint_left.x = mean([joint_center.x, node.x(node.id == joint.x_neg(i),:)]);
                        joint_left.y = joint_center.y;
                        joint_left.z = joint_center.z;
                        fprintf(fileID,'node %i %f %f %f \n',44000+i,joint_left.x,joint_left.y,joint_left.z);
                        ele_props = ele_props_table(ele_props_table.id == bm_left.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),44000+i,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+7,44000+i,joint_center_node_beams);
                    end
                    % Beam Right
                    if ~isempty(bm_right)
                        joint_right.x = mean([joint_center.x, node.x(node.id == joint.x_pos(i),:)]);
                        joint_right.y = joint_center.y;
                        joint_right.z = joint_center.z;
                        fprintf(fileID,'node %i %f %f %f \n',45000+i,joint_right.x,joint_right.y,joint_right.z);
                        ele_props = ele_props_table(ele_props_table.id == bm_right.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+2,45000+i,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+8,joint_center_node_beams,45000+i);
                    end
                    % Column Offsets
                    % Columns Below
                    joint_down.x = joint_center.x;
                    joint_down.y = mean([joint_center.y, node.y(node.id == joint.y_neg(i),:)]);
                    joint_down.z = joint_center.z;
                    fprintf(fileID,'node %i %f %f %f \n',54000+i,joint_down.x,joint_down.y,joint_down.z);
                    if ~isempty(col_low)
                        ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i);
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+9,54000+i,joint_center_node_columns);
                    % Columns Above
                    joint_up.x = joint_center.x;
                    joint_up.y = mean([joint_center.y, node.y(node.id == joint.y_pos(i),:)]);
                    joint_up.z = joint_center.z;
                    fprintf(fileID,'node %i %f %f %f \n',55000+i,joint_up.x,joint_up.y,joint_up.z);
                    if ~isempty(col_high)
                        ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i));
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+10,joint_center_node_columns,55000+i);
                end
                % Out of plane beam offsets
%                 fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*1000000+5,joint.z_neg(i),40000+i);
%                 fprintf(fileID,'element elasticBeamColumn %d %d %d 1000. 99999999. 99999999. 99999. 200000. 200000. 3 \n',joint.id(i)*1000000+6,40000+i,joint.z_pos(i));
                % Scissor Hinge
                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n', joint.id(i)+10000, joint_center_node_beams, joint_center_node_columns, joint.id(i)+10000); % Currently assumes x direction joint
                    fprintf(fileID,'equalDOF %i %i 1 2 3 4 5 \n', 40000+i, 50000+i);
                end
            elseif analysis.joint_model == 2  % Joint 3D
                % element Joint3D %tag %Nx- %Nx+ %Ny- %Ny+ %Nz- %Nz+ %Nc %MatX %MatY %MatZ %LrgDspTag
                fprintf(fileID,'element Joint3D %i %i %i %i %i %i %i %i %i 1 1 0 \n', 10000+i, joint.x_neg(i), joint.x_pos(i), joint.y_neg(i), joint.y_pos(i), joint.z_neg(i), joint.z_pos(i), 10000+i, 10000+i); 
            end
        end
    end
end

%% Define Rigid Slabs
if analysis.rigid_diaphram
    for s = 1:height(story)
%         if story.id(s) > 0
            slab_nodes_at_story = node.id(node.on_slab == 1 & node.story == story.id(s))';
            fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(slab_nodes_at_story));
%         end
    end
end
                
%% Define Plastic Hinges
if height(hinge) > 0
    % Load linear element table
    if analysis.nonlinear ~= 0 && analysis.model_type == 2 % Nonlinear MDOF
        if exist([read_dir_analysis filesep 'element_analysis.mat'],'file')
            element_analysis_temp = load([read_dir_analysis filesep 'element_analysis.mat']);
            element_analysis = element_analysis_temp.element;
        else
            [ element_analysis, joint ] = main_element_capacity( story, ele_props_table, element, analysis, joint, read_dir_analysis );
            [ element_analysis, ~ ] = main_hinge_properties( ele_props_table, element_analysis, joint );
        end
    else % SDOF
        element_analysis = element;
    end
    
    for i = 1:height(hinge)
        hin = hinge(i,:);
        ele_hinge_id = element.id(end) + hin.id; % Element ID associated with this hinge
        if strcmp(hin.type,'foundation')
            if strcmp(hin.direction,'pile')
                foundation_rot_stiff = 6920602859; % Force-Displacement rotational Stiffness of Bundle of Piles
                foundation_lat_stiff = 11026800; % Force-Displacement lateral Stiffness of Bundle of Piles
                foundation_axial_stiff = 5884866; % Force-Displacement lateral Stiffness of Bundle of Piles
            elseif strcmp(hin.direction,'wall')
                foundation_rot_stiff = 9999999999; % essentially rigid rotational stiffness
                foundation_lat_stiff = 1; % essentially no lateral stiffness
                foundation_axial_stiff = 153876; % axial stiffness of the wall with a gravity load sitting on the soil
            end
%             pile_lat_stiff = 1; % Force-Displacement lateral Stiffness of Bundle of Piles
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id,foundation_rot_stiff); % Elastic Material
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id+8000,foundation_lat_stiff); % Elastic Material
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id+9000,foundation_axial_stiff); % Elastic Material
            if strcmp(dimension,'3D') % Input as rotational for now
%                 fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 %i 2 %i -dir 1 2 3 4 5 6 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id, ele_hinge_id);
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i %i %i 2 %i -dir 1 2 3 4 5 6 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id+8000, ele_hinge_id+9000, ele_hinge_id+8000, ele_hinge_id, ele_hinge_id);
            else
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i %i -dir 1 2 3 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id+8000, ele_hinge_id+9000, ele_hinge_id);
            end
        else
            ele_side = num2str(hin.ele_side);
            ele = element(element.id == hin.element_id,:);
            if analysis.model_type == 2 % MDOF
                ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
            else % SDOF
                ele_props = element;
            end

            if analysis.nonlinear ~= 0
                hinge_props = element_analysis(element_analysis.id == hin.element_id,:);
            end
            if strcmp(ele_props.type,'beam') || strcmp(ele_props.type,'column')  % IMK Rotational Hinge
                % Define backbone curve
                if strcmp(hin.direction,'primary')
                    [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', hinge_props.(['Mn_pos_' ele_side]), hinge_props.(['Mn_neg_' ele_side]), hinge_props.(['Mp_pos_' ele_side]), hinge_props.(['Mp_neg_' ele_side]), ele.length, ele_props.e, ele_props.iz, hinge_props.(['a_hinge_' ele_side]), hinge_props.(['b_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), analysis.hinge_stiff_mod, 0.1, hinge_props.(['critical_mode_' ele_side]) );
                elseif strcmp(hin.direction,'oop')
                    [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', hinge_props.(['Mn_oop_' ele_side]), hinge_props.(['Mn_oop_' ele_side]), hinge_props.(['Mp_oop_' ele_side]), hinge_props.(['Mp_oop_' ele_side]), ele.length, ele_props.e, ele_props.iy, hinge_props.(['a_hinge_oop_' ele_side]), hinge_props.(['b_hinge_oop_' ele_side]), hinge_props.(['c_hinge_oop_' ele_side]), analysis.hinge_stiff_mod, 0.1, hinge_props.(['critical_mode_oop_' ele_side]) );
                end

                % Define IMK Parameters
                K0 = moment_vec_pos(1)/rot_vec_pos(1);
                as_sping_pos = (moment_vec_pos(2)-moment_vec_pos(1))/(rot_vec_pos(2)-rot_vec_pos(1))/K0;
                as_sping_neg = (moment_vec_neg(2)-moment_vec_neg(1))/(rot_vec_neg(2)-rot_vec_neg(1))/K0;
                theta_p_pos = rot_vec_pos(2)-rot_vec_pos(1);
                theta_p_neg = rot_vec_neg(2)-rot_vec_neg(1);
                residual_strength = 0.05; % fix to 5%
                q_ult_pos = moment_vec_pos(2)/moment_vec_pos(1);
                q_ult_neg = moment_vec_neg(2)/moment_vec_neg(1);
                theta_pc_pos = rot_vec_pos(3) - rot_vec_pos(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(rot_vec_pos(3) - rot_vec_pos(2))/(q_ult_pos-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
                theta_pc_neg = rot_vec_neg(3) - rot_vec_neg(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(rot_vec_neg(3) - rot_vec_neg(2))/(q_ult_neg-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
%                 if analysis.type == 1 % Dynamic
                    end_rot = 0.999999; % Keep residual strength forever
%                 else
%                     end_rot = rot_vec_pos(end);
%                 end
                % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%                 fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), theta_p_pos, theta_p_neg, theta_pc_pos, theta_pc_neg, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), end_rot, end_rot);

                % uniaxialMaterial IMKPeakOriented $Mat_Tag $Ke $Up_pos $Upc_pos $Uu_pos $Fy_pos $FmaxFy_pos $FresFy_pos $Up_neg $Upc_neg $Uu_neg $Fy_neg $FmaxFy_neg $FresFy_neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $D_pos $D_neg
                fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', ele_hinge_id, K0, theta_p_pos, theta_pc_pos, end_rot, moment_vec_pos(1), moment_vec_pos(2)/moment_vec_pos(1), residual_strength, theta_p_neg, theta_pc_neg, end_rot, moment_vec_neg(1), moment_vec_neg(2)/moment_vec_neg(1), residual_strength);
               
                % uniaxialMaterial Bilin $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg <$nFactor>
%                 fprintf(fileID,'uniaxialMaterial Bilin %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc_pos, theta_pc_neg, residual_strength, residual_strength, end_rot, end_rot);

                % Create Zero Length Elements and Equal DOF constraints
                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                if strcmp(ele.direction,'x') && strcmp(hin.direction,'oop') % Out of plane for the X direction  
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                elseif strcmp(ele.direction,'x') % In plane for the X direction 
                    if strcmp(ele_props.type,'column')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir %i \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id, rot_dof_x);
                        fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 -dir 1 2 3 5 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                    elseif strcmp(ele_props.type,'beam')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir %i \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id, rot_dof_x);
                        fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 2 -dir 1 2 3 4 5 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                    end
                elseif strcmp(ele.direction,'z') && strcmp(hin.direction,'oop') % Out of plane for the Z direction  
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                elseif strcmp(ele.direction,'z') % In plane for the Z direction 
                    if strcmp(ele_props.type,'column')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                        fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 -dir 1 2 3 5 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                    elseif strcmp(ele_props.type,'beam')
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                        fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 2 -dir 1 2 3 5 6 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                    end
                end
                
            elseif strcmp(ele_props.type,'wall') && ~analysis.fiber_walls % lumped plasticity wall models
                if analysis.nonlinear == 0 % Elastic Lateral Spring for shear deformations
                    elastic_shear_stiffness = ele_props.g*ele_props.av/ele.length;
                    % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
                    fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id,elastic_shear_stiffness); 
                elseif analysis.nonlinear == 1% Nonlinear
                    % Define backbone coordinates and IMK Hinges
                    if strcmp(hin.direction,'primary')
                        [ force_vec, disp_vec ] = fn_define_backbone_shear( hinge_props.(['Vn_' ele_side]), ele.length, ele_props.g, ele_props.av, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['d_hinge_' ele_side]), hinge_props.(['e_hinge_' ele_side]), hinge_props.(['f_hinge_' ele_side]), hinge_props.(['g_hinge_' ele_side])  );
                        K0 = force_vec(1)/disp_vec(1);
                        residual_strength = 0.05; % fix to 5% 
                        % Have it go past the shear kink with the initial stiffness and check how far it goes in post process
                        f_yield = force_vec(2);
                        f_ult_ratio = 1;
                        theta_p = disp_vec(3)-force_vec(2)/K0; % Correct theta P based initial elastic stiffness
                        theta_pc = disp_vec(4) - disp_vec(3) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(disp_vec(4) - disp_vec(3))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
%                         
                        % Have it go straight from yeild to residual
%                         f_yield = force_vec(1);
%                         f_ult_ratio = force_vec(2)/force_vec(1);
%                         theta_p = disp_vec(2)-disp_vec(1); % Theta P is the disp of the first kink
%                         theta_pc = disp_vec(4) - disp_vec(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(disp_vec(4) - disp_vec(2))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
%                         as_sping = (force_vec(2)-force_vec(1))/(disp_vec(2)-disp_vec(1))/K0;

                        % if analysis.type == 1 % Dynamic
                            end_disp = 999; % Keep residual strength forever
%                         else
%                             end_disp = disp_vec(end);
%                         end
                        % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%                         fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, K0, 0, 0, force_vec(2), -force_vec(2), theta_p, theta_p, theta_pc, theta_pc, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), end_disp, end_disp);
%                         fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping, as_sping, force_vec(1), -force_vec(1), theta_p, theta_p, theta_pc, theta_pc, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), end_disp, end_disp); % Keep residual strength forever
                        
                        % uniaxialMaterial IMKPeakOriented $Mat_Tag $Ke $Up_pos $Upc_pos $Uu_pos $Fy_pos $FmaxFy_pos $FresFy_pos $Up_neg $Upc_neg $Uu_neg $Fy_neg $FmaxFy_neg $FresFy_neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $D_pos $D_neg
                        fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', ele_hinge_id, K0, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength);
               
                      % uniaxialMaterial Bilin $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg <$nFactor>
%                         fprintf(fileID,'uniaxialMaterial Bilin %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping, as_sping, force_vec(1), -force_vec(1), theta_p, theta_p, theta_pc, theta_pc, residual_strength, residual_strength, end_disp, end_disp); % Keep residual strength forever
                    elseif strcmp(hin.direction,'oop')
                        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', hinge_props.(['Mn_oop_' ele_side]), hinge_props.(['Mn_oop_' ele_side]), hinge_props.(['Mp_oop_' ele_side]), hinge_props.(['Mp_oop_' ele_side]), ele.length, ele_props.e, ele_props.iy, hinge_props.(['a_hinge_oop_' ele_side]), hinge_props.(['b_hinge_oop_' ele_side]), hinge_props.(['c_hinge_oop_' ele_side]), analysis.hinge_stiff_mod, 0.1, hinge_props.(['critical_mode_oop_' ele_side]) );
                        K0 = moment_vec_pos(1)/rot_vec_pos(1);
                        residual_strength = 0.05; % fix to 5% 
                        as_sping_pos = (moment_vec_pos(2)-moment_vec_pos(1))/(rot_vec_pos(2)-rot_vec_pos(1))/K0;
                        as_sping_neg = (moment_vec_neg(2)-moment_vec_neg(1))/(rot_vec_neg(2)-rot_vec_neg(1))/K0;
                        theta_pc_pos = rot_vec_pos(3) - rot_vec_pos(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(rot_vec_pos(3) - rot_vec_pos(2))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
                        theta_pc_neg = rot_vec_neg(3) - rot_vec_neg(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(rot_vec_neg(3) - rot_vec_neg(2))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
%                         if analysis.type == 1 % Dynamic
                            end_rot = 0.999999;  % Keep residual strength forever
%                         else
%                             end_rot = rot_vec_pos(end);
%                         end
                        % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%                         fprintf(fileID,'uniaxialMaterial Bilin %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc_pos, theta_pc_neg, residual_strength, residual_strength, end_rot, end_rot); % Keep residual strength forever
%                         fprintf(fileID,'uniaxialMaterial ElasticPP %i %f %f \n',ele_hinge_id, Ko, rot_vec_pos(1));

                        % First story Hinges modified to match fiber model
                        if hin.story == 1
                            fprintf(fileID,'uniaxialMaterial ElasticBilin %i %f %f %f \n',ele_hinge_id+15000, K0, K0/1500, rot_vec_pos(1)*1.2);
                        else % second story Hinges modified to match fiber model
                            fprintf(fileID,'uniaxialMaterial ElasticBilin %i %f %f %f \n',ele_hinge_id+15000, K0, K0/3500, rot_vec_pos(1)*0.87);
                        end
                        % uniaxialMaterial MinMax $matTag $otherTag <-min $minStrain> <-max $maxStrain>
                        fprintf(fileID,'uniaxialMaterial MinMax %i %i -min %f -max %f \n', ele_hinge_id, ele_hinge_id+15000, -rot_vec_neg(3), rot_vec_pos(3));
                        
                        % Essentially Pinned OOP Wall
%                         fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id,1);
                    end       
                end
                
                % Create Zero Length Element (Currently does not work for 2D)
                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                if strcmp(ele.direction,'x') && strcmp(hin.direction,'oop') % Out of plane for the X direction  
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                elseif strcmp(ele.direction,'x') % In plane for the X direction (assume shear)
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 1 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                    fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 222222 222222 -dir 2 3 5 6 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                elseif strcmp(ele.direction,'z') && strcmp(hin.direction,'oop') % Out of plane for the Z direction  
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                elseif strcmp(ele.direction,'z') % In plane for the Z direction (assume shear)
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id);
                    fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 222222 222222 -dir 1 2 4 5 \n',50000 + ele_hinge_id, hin.node_1, hin.node_2);
                end          
            end
        end
    end
end

% Print model to file 
% fprintf(fileID,'print -file %s/model.txt \n',output_dir);

if ~analysis.suppress_outputs
    fprintf(fileID,'puts "Model Build Complete" \n');
end

% Close File
fclose(fileID);

end

