function [ joint_ele_ids ] = fn_define_model( write_dir, node, element, joint, hinge, analysis, dimension, story, read_dir_analysis, model, ele_props_table )
%UNTITLED6 Summary of this function goes here

%% Import Tools
import asce_41.*
import build_model.fn_node_exist

%% Load element properties table
% if analysis.model_type ~= 3
% %     ele_props_table = readtable([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'],'Sheet','element'); % for archetype models, the model properties are already in the table
if analysis.model_type ~= 3
    ele_props_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
end

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
        else %MDOF or archetype
            if analysis.nonlinear ~= 0 && ~element.elastic(i) % Nonlinear Element
                if strcmp(analysis.nonlinear_type,'lumped')
                    % Lumped Plasticity Model of beams and columns
                    
                    % Add stiffness to element to account for two springs, from appendix B of Ibarra and Krawinkler 2005
                    if analysis.model_type == 3 % Archetype model
                        Iz_ele = ele_props.iz_model*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod); 
                    else
                        Iz_ele = ele_props.iz*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod);
                    end

                    if strcmp(dimension,'2D')
                        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,Iz_ele,geotransf);
                        % element ElasticTimoshenkoBeam $eleTag $iNode $jNode $E $G $A $Iz $Avy $transfTag <-mass $massDens> <-cMass>
    %                     fprintf(fileID,'element ElasticTimoshenkoBeam %i %i %i %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.e,ele_props.g,ele_props.a,Iz_ele,(5/6)*ele_props.a,geotransf);
                    elseif strcmp(dimension,'3D')
                        Iy_ele = ele_props.iy*((analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod);

                        % element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,Iy_ele,Iz_ele,geotransf);
                    end
                    
                elseif strcmp(analysis.nonlinear_type,'fiber')
                    % Fiber model of beams and columns
                    
                    % Fixed steel parameters
                    fy_e = 70000;
                    fy_h = fy_e;
                    Es = 29000000;
                    e_su = 0.09; % From ATC 114 vol 3 apendix B
                    
                    % Fixed concrete parameters
                    fp_co = ele_props.fc_e;
                    ec_o = 0.002;
                    ec_sp = 2*ec_o;
                    ec_t = 0.0001;
                    fc_t = 7.5*sqrt(ele_props.fc_n);
                    
                    % Caluclate concrete parameters from ATC 114
                    ke = (1-sum(ele_props.w_p)/(6*ele_props.b_c*ele_props.d_c))*(1-ele_props.s_p/(2*ele_props.b_c))*(1-ele_props.s_p/(2*ele_props.d_c))/(1-ele_props.rho_cc); % equation 5-9 from ATC 114 vol 3
                    fp_l = 0.5*ke*ele_props.rho_s*fy_h; % equation 5-8 from ATC 114 vol 3
                    fp_cc = fp_co*(-1.254 + 2.254*sqrt(1 + 7.94*fp_l/fp_co) - 2*fp_l/fp_co); % equation 5-10 from ATC 114 vol 3
                    ec_c = ec_o*(1 + 5*(fp_cc/fp_co - 1)); % equation 5-7 from ATC 114 vol 3
                    ec_u = 0.004 + 1.4*ele_props.rho_s*fy_h*e_su/fp_cc; % equation 5-11 from ATC 114 vol 3
                    
                    % uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
                    fprintf(fileID,'uniaxialMaterial Steel02 %i %f %f 0.01 15.0 0.925 0.15 \n', 1000 + element.id(i), fy_e, Es);

                    % uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
                    fprintf(fileID,'uniaxialMaterial Concrete04 %i %f %f %f %f %f %f \n', 2000 + element.id(i), -fp_co, -ec_o, -ec_sp, ele_props.e, fc_t, ec_t); % UnConfined concrete

                    % uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
                    fprintf(fileID,'uniaxialMaterial Concrete04 %i %f %f %f %f %f %f \n', 3000 + element.id(i), -fp_cc, -ec_c, -ec_u, ele_props.e, fc_t, ec_t); % Confined concrete
                        
                    % section Fiber $secTag <-GJ $GJ> {
                    fprintf(fileID,'section Fiber %i { \n',element.id(i));
                    
                    As = str2double(strsplit(strrep(strrep(ele_props.As{1},'[',''),']','')));
                    num_bars = str2double(strsplit(strrep(strrep(ele_props.n_b{1},'[',''),']','')));
                    depth_bars = str2double(strsplit(strrep(strrep(ele_props.As_d{1},'[',''),']','')));
    
                   
                    % Web/Core
                    % patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), 1, 1, ele_props.d_c/2, -ele_props.w/2, ele_props.h/2, ele_props.w/2); % top unconfined
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), 1, 1, -ele_props.h/2, -ele_props.w/2, -ele_props.d_c/2, ele_props.w/2); % bottom unconfined
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), floor(ele_props.d_c), 1, -ele_props.d_c/2, -ele_props.w/2,  ele_props.d_c/2, -ele_props.b_c/2); % left unconfined
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), floor(ele_props.d_c), 1, -ele_props.d_c/2, ele_props.b_c/2, ele_props.d_c/2, ele_props.w/2); % right unconfined
                    fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 3000 + element.id(i), floor(ele_props.d_c), 1, -ele_props.d_c/2, -ele_props.b_c/2,  ele_props.d_c/2, ele_props.b_c/2); % confined core
                    
%                     fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 3000 + element.id(i), ele_props.h, 1, -ele_props.h/2, -ele_props.w/2, ele_props.h/2, ele_props.w/2); % confined core
                        
                    if strcmp(element.type{i},'beam') && ele_props.d_flange > 0
                        % Flange to the left and right
                        % patch rect $matTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ
                        fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), 1, 1, ele_props.h/2-ele_props.d_flange, -ele_props.b_eff/2, ele_props.h/2, -ele_props.w/2);
                        fprintf(fileID,'patch rect %i %i %i %f %f %f %f \n', 2000 + element.id(i), 1, 1, ele_props.h/2-ele_props.d_flange, ele_props.w/2, ele_props.h/2, ele_props.b_eff/2 );
                    end

                    % layer straight $matTag $numFiber $areaFiber $yStart $zStart $yEnd $zEnd
                    yStart = -ele_props.w/2 + ele_props.clear_cover;
                    yEnd = ele_props.w/2 - ele_props.clear_cover;
                    for row = 1:length(As)
                        if As(row) > 0
                            z_coor = ele_props.h/2 - depth_bars(row);
                            fprintf(fileID,'layer straight %i %i %f %f %f %f %f \n', 1000 + element.id(i), num_bars(row), As(row), yStart, z_coor, yEnd, z_coor);
                        end
                    end

                    fprintf(fileID,'} \n');

                    % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>
                    fprintf(fileID,'element forceBeamColumn %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),5,element.id(i),geotransf);  
                else
                    error('Nonlinear Model Type Not Recognized')
                end
            else
                if strcmp(dimension,'2D')
                    % element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
                elseif strcmp(dimension,'3D')
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
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.iz,geotransf);
                elseif strcmp(dimension,'3D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,geotransf);
                end
            elseif analysis.fiber_walls % model with fiber elements
                % uniaxialMaterial Steel02 $matTag $Fy $E $b $R0 $cR1 $cR2 <$a1 $a2 $a3 $a4 $sigInit>
                fprintf(fileID,'uniaxialMaterial Steel02 %i %f %f 0.01 15.0 0.925 0.15 \n', 1000 + element.id(i), ele_props.fy_e, ele_props.Es);
                % uniaxialMaterial Concrete04 $matTag $fc $ec $ecu $Ec <$fct $et> <$beta>
                fprintf(fileID,'uniaxialMaterial Concrete04 %i %f -0.002 -0.006 %f %f 0.0001 \n',2000 + element.id(i), -ele_props.fc_e, ele_props.e, 7.5*sqrt(ele_props.fc_n));                
                
                if analysis.fiber_walls == 1 % center-line beam column fiber model
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

                    fprintf(fileID,'} \n');

                    % element forceBeamColumn $eleTag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens> <-iter $maxIters $tol> <-integration $intType>
                    fprintf(fileID,'element forceBeamColumn %i %i %i %i %i %i \n',element.id(i),element.node_1(i),element.node_2(i),5,element.id(i),geotransf);  
                elseif analysis.fiber_walls == 2 % MVLEM
                    [ force_vec, disp_vec ] = fn_define_backbone_shear( 1.34e6, element.length(i), ele_props.g, ele_props.av, 0.2, 1, 2, 0.6, 0.4 );
                    K0 = force_vec(1)/disp_vec(1);
                    residual_strength = 0.05; % fix to 5% 

                    % Have it go straight from yeild to residual
                    f_yield = force_vec(1);
                    f_ult_ratio = force_vec(2)/force_vec(1);
                    theta_p = disp_vec(2)-disp_vec(1); % Theta P is the disp of the first kink
                    theta_pc = disp_vec(4) - disp_vec(2) + (0.2-residual_strength)*(disp_vec(4) - disp_vec(2))/(1-0.2); % theta pc defined all the way to zero where b defined to residual kink

                    end_disp = 999; % Keep residual strength forever

                    % uniaxialMaterial IMKPeakOriented $Mat_Tag $Ke $Up_pos $Upc_pos $Uu_pos $Fy_pos $FmaxFy_pos $FresFy_pos $Up_neg $Upc_neg $Uu_neg $Fy_neg $FmaxFy_neg $FresFy_neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $D_pos $D_neg
                    fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 100.0 100.0 100.0 100.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', 3000 + element.id(i), K0, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength);

                    
                    % element MVLEM $eleTag $Dens $iNode $jNode $m $c -thick {Thicknesses} -width {Widths} -rho {Reinforcing_ratios} -matConcrete {Concrete_tags} -matSteel {Steel_tags} -matShear {Shear_tag}
                    fprintf(fileID,'element MVLEM %i 0.0 %i %i 6 0.4 -thick 12.0 12.0 12.0 12.0 12.0 12.0 -width 37.5 37.5 37.5 37.5 37.5 37.5 -rho 0.002 0.002 0.002 0.002 0.002 0.002 -matConcrete %i %i %i %i %i %i -matSteel %i %i %i %i %i %i -matShear %i \n',element.id(i),element.node_1(i),element.node_2(i), 2000 + element.id(i), 2000 + element.id(i), 2000 + element.id(i), 2000 + element.id(i), 2000 + element.id(i), 2000 + element.id(i), 1000 + element.id(i), 1000 + element.id(i), 1000 + element.id(i), 1000 + element.id(i), 1000 + element.id(i), 1000 + element.id(i), 3000 + element.id(i));
%                     fprintf(fileID,'Element MVLEM 1 0.0 1 2 8 0.4 -thick 4 4 4 4 4 4 4 4 -width 7.5 1.5 7.5 7.5 7.5 7.5 1.5 7.5 -rho 0.0293 0.0 0.0033 0.0033 0.0033 0.0033 0.0 0.0293 -matConcrete 3 4 4 4 4 4 4 3 -matSteel 1 2 2 2 2 2 2 1 -matShear 5 \n')
                end
            end
        end

    %% Truss Assigment
    elseif strcmp(element.type{i},'truss')
        % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
        fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',element.id(i),ele_props.e);
        % element truss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
        fprintf(fileID,'element truss %i %i %i %f %i -rho 0.0 -cMass 0 -doRayleigh 0 \n',element.id(i),element.node_1(i),element.node_2(i),ele_props.a,element.id(i));
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
        if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1 && joint.story(i) <= analysis.stories_nonlinear && joint.story(i) > analysis.stories_nonlinear_low % Nonlinear Joints
            [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', joint_analysis.Mn(i), joint_analysis.Mn(i), joint_analysis.Mn(i), joint_analysis.Mn(i), joint_analysis.h(i), joint_analysis.e(i), joint_analysis.iz(i), joint_analysis.a_hinge(i), joint_analysis.b_hinge(i), joint_analysis.c_hinge(i), 100, 0, 'shear' );
            K0 = moment_vec_pos(1)/rot_vec_pos(1);
            theta_p_pos = rot_vec_pos(2)-rot_vec_pos(1);
            theta_p_neg = rot_vec_neg(2)-rot_vec_neg(1);
            residual_strength = 0.05; % fix to 5%
            q_ult_pos = moment_vec_pos(2)/moment_vec_pos(1);
            q_ult_neg = moment_vec_neg(2)/moment_vec_neg(1);
            theta_pc_pos = rot_vec_pos(3) - rot_vec_pos(2) + (joint_analysis.c_hinge(i)-residual_strength)*(rot_vec_pos(3) - rot_vec_pos(2))/(q_ult_pos-joint_analysis.c_hinge(i)); % theta pc defined all the way to zero where b defined to residual kink
            theta_pc_neg = rot_vec_neg(3) - rot_vec_neg(2) + (joint_analysis.c_hinge(i)-residual_strength)*(rot_vec_neg(3) - rot_vec_neg(2))/(q_ult_neg-joint_analysis.c_hinge(i)); % theta pc defined all the way to zero where b defined to residual kink
            end_rot = 0.999999; % Keep residual strength forever

            % uniaxialMaterial IMKPeakOriented $Mat_Tag $Ke $Up_pos $Upc_pos $Uu_pos $Fy_pos $FmaxFy_pos $FresFy_pos $Up_neg $Upc_neg $Uu_neg $Fy_neg $FmaxFy_neg $FresFy_neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $D_pos $D_neg
            fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 100.0 100.0 100.0 100.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', joint.id(i)+10000, K0, theta_p_pos, theta_pc_pos, end_rot, moment_vec_pos(1), moment_vec_pos(2)/moment_vec_pos(1), residual_strength, theta_p_neg, theta_pc_neg, end_rot, moment_vec_neg(1), moment_vec_neg(2)/moment_vec_neg(1), residual_strength);
        end
        if analysis.joint_model == 1 % Elastic beam column elements
            joint_ele_ids = [(1:6)+joint.id(i)*1000000,joint_ele_ids];
            joint_center.x = node.x(node.id == joint.y_pos(i),:);
            joint_center.y = node.y(node.id == joint.x_pos(i),:);
            joint_center.z = node.z(node.id == joint.x_pos(i),:);
            
            exist_node_filt = node.x == joint_center.x & node.y == joint_center.y & node.z == joint_center.z; 
            if sum(exist_node_filt) == 0 % no existing nodes, so create a new one(s)
                if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1  && joint.story(i) <= analysis.stories_nonlinear && joint.story(i) > analysis.stories_nonlinear_low %  Explicit Joint Model
                    joint_center_node_beams = 40000+i;
                    joint_center_node_columns = 50000+i;
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'node %i %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y);
                        fprintf(fileID,'node %i %f %f \n',joint_center_node_columns,joint_center.x,joint_center.y);
                    else
                        fprintf(fileID,'node %i %f %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y,joint_center.z);
                        fprintf(fileID,'node %i %f %f %f \n',joint_center_node_columns,joint_center.x,joint_center.y,joint_center.z);
                    end
                else
                    joint_center_node_beams = 40000+i;
                    joint_center_node_columns = 40000+i;
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'node %i %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y);
                    else
                        fprintf(fileID,'node %i %f %f %f \n',joint_center_node_beams,joint_center.x,joint_center.y,joint_center.z);
                    end
                end
            elseif sum(exist_node_filt) == 1 % existing node(s), use those
                % Use existing node
                if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1  && joint.story(i) <= analysis.stories_nonlinear && joint.story(i) > analysis.stories_nonlinear_low %  Explicit Joint Model
                    joint_center_node_beams = node.id(exist_node_filt);
                    % defin new node for zero length element to connet
                    joint_center_node_columns = 50000 + node.id(exist_node_filt);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'node %i %f %f \n',joint_center_node_columns,joint_center.x,joint_center.y);
                    else
                        fprintf(fileID,'node %i %f %f %f \n',joint_center_node_columns,joint_center.x,joint_center.y,joint_center.z);
                    end
                else
                    joint_center_node_beams = node.id(exist_node_filt);
                    joint_center_node_columns = node.id(exist_node_filt);
                end
            else
                % found multiple nodes! Should not get here
                error('Multiple Nodes Found. Revise code')
            end

            % Joint Fixity Condition
            if exist('joint_analysis','var')
                joint_rigidity = joint_analysis.implicit_stiff(i);
            else
                % Default assumption to Rigid columns sections in joints
                joint_rigidity = 1;
            end

            % Grab element info
            bm_left = element(element.id == joint.beam_left(i),:);
            bm_right = element(element.id == joint.beam_right(i),:);
            col_low = element(element.id == joint.column_low(i),:);
            col_high = element(element.id == joint.column_high(i),:);

            % Primary Beam Offsets
            if joint_rigidity == 1 % Column Offsets Rigid
                % Beam Offsets
                if ~isempty(bm_left)
                    ele_props = ele_props_table(ele_props_table.id == bm_left.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams,ele_props.a,ele_props.e,ele_props.iz,2);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                    end
                end
                if ~isempty(bm_right)
                    ele_props = ele_props_table(ele_props_table.id == bm_right.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.iz,2);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                    end
                end
                % Column Offsets
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                end
            elseif joint_rigidity == 2 % Beam Offsets Rigid
                % Primary Beam Offsets
                if ~isempty(bm_left)
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node_beams);
                    end
                end
                if ~isempty(bm_right)
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i));
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node_beams,joint.x_pos(i));
                    end
                end
                % Column Offsets
                if ~isempty(col_low)
                    ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns,ele_props.a,ele_props.e,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    end
                else
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node_columns);
                    end
                end
                if ~isempty(col_high)
                    ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    end
                else
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node_columns,joint.y_pos(i));
                    end
                end
            elseif joint_rigidity == 3 % Both Offsets are Half Rigid
                % Primary Beam Offsets
                % Beam Left
                if ~isempty(bm_left)
                    joint_left.x = mean([joint_center.x, node.x(node.id == joint.x_neg(i),:)]);
                    joint_left.y = joint_center.y;
                    joint_left.z = joint_center.z;
                    ele_props = ele_props_table(ele_props_table.id == bm_left.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'node %i %f %f \n',44000+i,joint_left.x,joint_left.y);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),44000+i,ele_props.a,ele_props.e,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+7,44000+i,joint_center_node_beams);
                    else
                        fprintf(fileID,'node %i %f %f %f \n',44000+i,joint_left.x,joint_left.y,joint_left.z);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),44000+i,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+7,44000+i,joint_center_node_beams);
                    end
                end
                % Beam Right
                if ~isempty(bm_right)
                    joint_right.x = mean([joint_center.x, node.x(node.id == joint.x_pos(i),:)]);
                    joint_right.y = joint_center.y;
                    joint_right.z = joint_center.z;
                    ele_props = ele_props_table(ele_props_table.id == bm_right.ele_id,:);
                    if strcmp(dimension,'2D')
                        fprintf(fileID,'node %i %f %f \n',45000+i,joint_right.x,joint_right.y);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+2,45000+i,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+8,joint_center_node_beams,45000+i);
                    else
                        fprintf(fileID,'node %i %f %f %f \n',45000+i,joint_right.x,joint_right.y,joint_right.z);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+2,45000+i,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+8,joint_center_node_beams,45000+i);
                    end
                end
                % Column Offsets
                % Columns Below
                joint_down.x = joint_center.x;
                joint_down.y = mean([joint_center.y, node.y(node.id == joint.y_neg(i),:)]);
                joint_down.z = joint_center.z;
                if strcmp(dimension,'2D')
                    fprintf(fileID,'node %i %f %f \n',54000+i,joint_down.x,joint_down.y);
                    if ~isempty(col_low)
                        ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i,ele_props.a,ele_props.e,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i);
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+9,54000+i,joint_center_node_columns);
                else
                     fprintf(fileID,'node %i %f %f %f \n',54000+i,joint_down.x,joint_down.y,joint_down.z);
                    if ~isempty(col_low)
                        ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),54000+i);
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+9,54000+i,joint_center_node_columns);
                end
                % Columns Above
                joint_up.x = joint_center.x;
                joint_up.y = mean([joint_center.y, node.y(node.id == joint.y_pos(i),:)]);
                joint_up.z = joint_center.z;
                if strcmp(dimension,'2D')
                    fprintf(fileID,'node %i %f %f \n',55000+i,joint_up.x,joint_up.y);
                    if ~isempty(col_high)
                        ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i));
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+10,joint_center_node_columns,55000+i);
                else
                    fprintf(fileID,'node %i %f %f %f \n',55000+i,joint_up.x,joint_up.y,joint_up.z);
                    if ~isempty(col_high)
                        ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                        fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                    else
                        fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,55000+i,joint.y_pos(i));
                    end
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+10,joint_center_node_columns,55000+i);
                end
            end
            % Scissor Hinge
            %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
            if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1 && joint.story(i) <= analysis.stories_nonlinear && joint.story(i) > analysis.stories_nonlinear_low
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 %i -dir 1 2 3 \n', joint.id(i)+10000, joint_center_node_beams, joint_center_node_columns, joint.id(i)+10000); % Currently assumes x direction joint
                else
                    fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n', joint.id(i)+10000, joint_center_node_beams, joint_center_node_columns, joint.id(i)+10000); % Currently assumes x direction joint
                    fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 2 -dir 1 2 3 4 5 \n',joint.id(i)+20000, joint_center_node_beams, joint_center_node_columns);
                end
            end
        elseif analysis.joint_model == 2  % Joint 2D and 3D
            if strcmp(dimension,'2D')
                % element Joint2D $eleTag $Nd1 $Nd2 $Nd3 $Nd4 $NdC <$Mat1 $Mat2 $Mat3 $Mat4> $MatC $LrgDspTag
                fprintf(fileID,'element Joint2D %i %i %i %i %i %i 1 0 \n', 10000+i, joint.y_pos(i), joint.x_neg(i), joint.y_neg(i), joint.x_pos(i), 10000+i);
            else
                % element Joint3D %tag %Nx- %Nx+ %Ny- %Ny+ %Nz- %Nz+ %Nc %MatX %MatY %MatZ %LrgDspTag
                fprintf(fileID,'element Joint3D %i %i %i %i %i %i %i %i %i 1 1 0 \n', 10000+i, joint.x_neg(i), joint.x_pos(i), joint.y_neg(i), joint.y_pos(i), joint.z_neg(i), joint.z_pos(i), 10000+i, 10000+i); 
            end
        elseif analysis.joint_model == 0  % Centerline Model
            joint_ele_ids = [(1:6)+joint.id(i)*1000000,joint_ele_ids];
            joint_center.x = node.x(node.id == joint.y_pos(i),:);
            joint_center.y = node.y(node.id == joint.x_pos(i),:);
            joint_center.z = node.z(node.id == joint.x_pos(i),:);
            exist_node_filt = node.x == joint_center.x & node.y == joint_center.y & node.z == joint_center.z; 
            if sum(exist_node_filt) == 0 
                % No current node, therefore create a new one
                joint_center_node = 40000+i;
                if strcmp(dimension,'2D')
                    fprintf(fileID,'node %i %f %f \n',joint_center_node,joint_center.x,joint_center.y);
                else
                    fprintf(fileID,'node %i %f %f %f \n',joint_center_node,joint_center.x,joint_center.y,joint_center.z);
                end
            elseif sum(exist_node_filt) == 1 
                % Use existing node
                joint_center_node = node.id(exist_node_filt);
            else
                % found multiple nodes! Should not get here
                error('Multiple Nodes Found. Revise code')
            end
            
            

            bm_left = element(element.id == joint.beam_left(i),:);
            bm_right = element(element.id == joint.beam_right(i),:);
            col_low = element(element.id == joint.column_low(i),:);
            col_high = element(element.id == joint.column_high(i),:);
            if ~isempty(bm_left)
                ele_props = ele_props_table(ele_props_table.id == bm_left.ele_id,:);
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node,ele_props.a,ele_props.e,ele_props.iz,2);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                end
            end
            if ~isempty(bm_right)
                ele_props = ele_props_table(ele_props_table.id == bm_right.ele_id,:);
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.iz,2);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,2);
                end
            end
            if ~isempty(col_low)
                ele_props = ele_props_table(ele_props_table.id == col_low.ele_id,:);
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node,ele_props.a,ele_props.e,ele_props.iz,1);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node,ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                end
            end
            if ~isempty(col_high)
                ele_props = ele_props_table(ele_props_table.id == col_high.ele_id,:);
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %i \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.iz,1);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i %f %f %f %f %f %f %i \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i),ele_props.a,ele_props.e,ele_props.g,ele_props.j,ele_props.iy,ele_props.iz,1);
                end
            else
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                end
            end
        elseif analysis.joint_model == 3  % Rigid joint model
            joint_ele_ids = [(1:6)+joint.id(i)*1000000,joint_ele_ids];
            joint_center.x = node.x(node.id == joint.y_pos(i),:);
            joint_center.y = node.y(node.id == joint.x_pos(i),:);
            joint_center.z = node.z(node.id == joint.x_pos(i),:);
            
            exist_node_filt = node.x == joint_center.x & node.y == joint_center.y & node.z == joint_center.z; 
            if sum(exist_node_filt) == 0 
                % No current node, therefore create a new one
                joint_center_node = 40000+i;
                if strcmp(dimension,'2D')
                    fprintf(fileID,'node %i %f %f \n',joint_center_node,joint_center.x,joint_center.y);
                else
                    fprintf(fileID,'node %i %f %f %f \n',joint_center_node,joint_center.x,joint_center.y,joint_center.z);
                end
            elseif sum(exist_node_filt) == 1 
                % Use existing node
                joint_center_node = node.id(exist_node_filt);
            else
                % found multiple nodes! Should not get here
                error('Multiple Nodes Found. Revise code')
            end

            bm_left = element(element.id == joint.beam_left(i),:);
            bm_right = element(element.id == joint.beam_right(i),:);
            col_low = element(element.id == joint.column_low(i),:);
            col_high = element(element.id == joint.column_high(i),:);
            if ~isempty(bm_left)
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node);
                end
            end
            if ~isempty(bm_right)
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i));
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i));
                end
            end
            if ~isempty(col_low)
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node);
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node);
                end
            end
            if ~isempty(col_high)
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                end
            else
                if strcmp(dimension,'2D')
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                else
                    fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
                end
            end
            
            
%             joint_ele_ids = [(1:6)+joint.id(i)*1000000,joint_ele_ids];
%             joint_center.x = node.x(node.id == joint.y_pos(i),:);
%             joint_center.y = node.y(node.id == joint.x_pos(i),:);
%             joint_center.z = node.z(node.id == joint.x_pos(i),:);
%             joint_center_node = 40000+i;
%             if strcmp(dimension,'2D')
%                 fprintf(fileID,'node %i %f %f \n',joint_center_node,joint_center.x,joint_center.y);
%             else
%                 fprintf(fileID,'node %i %f %f %f \n',joint_center_node,joint_center.x,joint_center.y,joint_center.z);
%             end
%                 
%             bm_left = element(element.id == joint.beam_left(i),:);
%             bm_right = element(element.id == joint.beam_right(i),:);
%             col_low = element(element.id == joint.column_low(i),:);
%             col_high = element(element.id == joint.column_high(i),:);
%             
%             % Primary Beam Offsets
%             if ~isempty(bm_left)
%                 if strcmp(dimension,'2D')
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node);
%                 else
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+1,joint.x_neg(i),joint_center_node);
%                 end
%             end
%             if ~isempty(bm_right)
%                 if strcmp(dimension,'2D')
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i));
%                 else
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 2 \n',joint.id(i)*1000000+2,joint_center_node,joint.x_pos(i));
%                 end
%             end
%             
%             % Column Offsets
%             if ~isempty(col_low)
%                 if strcmp(dimension,'2D')
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node);
%                 else
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+3,joint.y_neg(i),joint_center_node);
%                 end
%             end
%             if ~isempty(col_high)
%                 if strcmp(dimension,'2D')
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
%                 else
%                     fprintf(fileID,'element elasticBeamColumn %i %i %i 1000. 99999999. 99999999. 99999. 200000. 200000. 1 \n',joint.id(i)*1000000+4,joint_center_node,joint.y_pos(i));
%                 end
%             end
        end
    end
end

%% Define Rigid Slabs
if strcmp(dimension,'3D') && analysis.rigid_diaphram
    for s = 1:height(story)
        slab_nodes_at_story = node.id(node.on_slab == 1 & node.story == story.id(s))';
        fprintf(fileID,'rigidDiaphragm 2 %s \n',num2str(slab_nodes_at_story));
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
    else % SDOF or Archetype
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
%                 foundation_axial_stiff = 153876; % axial stiffness of the wall with a gravity load sitting on the soil
                foundation_axial_stiff = 5884866; % Force-Displacement lateral Stiffness of Bundle of Piles
            end
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id,foundation_rot_stiff); % Elastic Material
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id+8000,foundation_lat_stiff); % Elastic Material
            fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id+9000,foundation_axial_stiff); % Elastic Material
            if strcmp(dimension,'3D') % Input as rotational for now
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i %i %i 2 %i -dir 1 2 3 4 5 6 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id+8000, ele_hinge_id+9000, ele_hinge_id+8000, ele_hinge_id, ele_hinge_id);
            else
                fprintf(fileID,'element zeroLength %i %i %i -mat %i %i %i -dir 1 2 3 \n',ele_hinge_id, hin.node_1, hin.node_2, ele_hinge_id+8000, ele_hinge_id+9000, ele_hinge_id);
            end
        else
            ele_side = num2str(hin.ele_side);
            ele = element(element.id == hin.element_id,:);
            if strcmp(ele_side,'1')
                nd_1 = hin.node_1;
                nd_2 = hin.node_2;
            elseif strcmp(ele_side,'2') % for side two, switch pos and negative to make sign correct
                nd_1 = hin.node_2;
                nd_2 = hin.node_1;    
            end

            if analysis.model_type == 2 % MDOF
                ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
            elseif analysis.model_type == 1 % SDOF
                ele_props = element;
            elseif analysis.model_type == 3 % archetype
                ele_props = element_analysis(element_analysis.id == hin.element_id,:);
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
                fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 100.0 100.0 100.0 100.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', ele_hinge_id, K0, theta_p_pos, theta_pc_pos, end_rot, moment_vec_pos(1), moment_vec_pos(2)/moment_vec_pos(1), residual_strength, theta_p_neg, theta_pc_neg, end_rot, moment_vec_neg(1), moment_vec_neg(2)/moment_vec_neg(1), residual_strength);

                % uniaxialMaterial Bilin $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg <$nFactor>
%                 fprintf(fileID,'uniaxialMaterial Bilin %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping_pos, as_sping_neg, moment_vec_pos(1), -moment_vec_neg(1), rot_vec_pos(2)-rot_vec_pos(1), rot_vec_neg(2)-rot_vec_neg(1), theta_pc_pos, theta_pc_neg, residual_strength, residual_strength, end_rot, end_rot);

                % Create Zero Length Elements and Equal DOF constraints
                %element zeroLength $eleTag $iNode $jNode -mat $matTag1 $matTag2 ... -dir $dir1 $dir2
                if strcmp(dimension,'2D')
                    if strcmp(ele.direction,'x') && strcmp(hin.direction,'primary') % Out of plane for the X direction  
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 3 \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id);
                        fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 -dir 1 2 \n',50000 + ele_hinge_id, nd_1, nd_2);
                    end
                elseif strcmp(dimension,'3D')
                    if strcmp(ele.direction,'x') && strcmp(hin.direction,'oop') % Out of plane for the X direction  
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id);
                    elseif strcmp(ele.direction,'x') % In plane for the X direction 
                        if strcmp(ele_props.type,'column')
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir %i \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id, rot_dof_x);
                            fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 -dir 1 2 3 5 \n',50000 + ele_hinge_id, nd_1, nd_2);
                        elseif strcmp(ele_props.type,'beam')
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir %i \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id, rot_dof_x);
                            fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 2 -dir 1 2 3 4 5 \n',50000 + ele_hinge_id, nd_1, nd_2);
                        end
                    elseif strcmp(ele.direction,'z') && strcmp(hin.direction,'oop') % Out of plane for the Z direction  
                        fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 6 \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id);
                    elseif strcmp(ele.direction,'z') % In plane for the Z direction 
                        if strcmp(ele_props.type,'column')
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id);
                            fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 -dir 1 2 3 5 \n',50000 + ele_hinge_id, nd_1, nd_2);
                        elseif strcmp(ele_props.type,'beam')
                            fprintf(fileID,'element zeroLength %i %i %i -mat %i -dir 4 \n',ele_hinge_id, nd_1, nd_2, ele_hinge_id);
                            fprintf(fileID,'element zeroLength %i %i %i -mat 1 1 1 2 2 -dir 1 2 3 5 6 \n',50000 + ele_hinge_id, nd_1, nd_2);
                        end
                    end
                end
                
            elseif strcmp(ele_props.type,'wall') && ~analysis.fiber_walls % lumped plasticity wall models
                if analysis.nonlinear == 0 % Elastic Lateral Spring for shear deformations
                    elastic_shear_stiffness = ele_props.g*ele_props.av/ele.length;
                    % uniaxialMaterial Elastic $matTag $E <$eta> <$Eneg>
                    fprintf(fileID,'uniaxialMaterial Elastic %i %f \n',ele_hinge_id,elastic_shear_stiffness); 
                    % define oop rotation rigidity that is not handled by
                    % the logic below for linear models
                    fprintf(fileID,'element zeroLength %i %i %i -mat 222222 -dir 6 \n',ele_hinge_id + 8675309, hin.node_1, hin.node_2);
                elseif analysis.nonlinear == 1% Nonlinear
                    % Define backbone coordinates and IMK Hinges
                    if strcmp(hin.direction,'primary')
                        [ force_vec, disp_vec ] = fn_define_backbone_shear( hinge_props.(['Vn_' ele_side]), ele.length, ele_props.g, ele_props.av, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['d_hinge_' ele_side]), hinge_props.(['e_hinge_' ele_side]), hinge_props.(['f_hinge_' ele_side]), hinge_props.(['g_hinge_' ele_side])  );
                        K0 = force_vec(1)/disp_vec(1);
                        residual_strength = 0.05; % fix to 5% 
%                         % Have it go past the shear kink with the initial stiffness for collapse assessment
%                         f_yield = force_vec(2);
%                         f_ult_ratio = 1;
%                         theta_p = disp_vec(3)-force_vec(2)/K0; % Correct theta P based initial elastic stiffness
%                         theta_pc = disp_vec(4) - disp_vec(3) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(disp_vec(4) - disp_vec(3))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink
%                          
                        % Have it go straight from yeild to residual
                        f_yield = force_vec(1);
                        f_ult_ratio = force_vec(2)/force_vec(1);
                        theta_p = disp_vec(2)-disp_vec(1); % Theta P is the disp of the first kink
                        theta_pc = disp_vec(4) - disp_vec(2) + (hinge_props.(['c_hinge_' ele_side])-residual_strength)*(disp_vec(4) - disp_vec(2))/(1-hinge_props.(['c_hinge_' ele_side])); % theta pc defined all the way to zero where b defined to residual kink

                        % if analysis.type == 1 % Dynamic
                            end_disp = 999; % Keep residual strength forever
%                         else
%                             end_disp = disp_vec(end);
%                         end
                        % uniaxialMaterial ModIMKPeakOriented $matTag $K0 $as_Plus $as_Neg $My_Plus $My_Neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
%                         fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, K0, 0, 0, force_vec(2), -force_vec(2), theta_p, theta_p, theta_pc, theta_pc, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), end_disp, end_disp);
%                         fprintf(fileID,'uniaxialMaterial ModIMKPeakOriented %i %f %f %f %f %f 10.0 10.0 10.0 10.0 1.0 1.0 1.0 1.0 %f %f %f %f %f %f %f %f 1.0 1.0 \n',ele_hinge_id, Ko, as_sping, as_sping, force_vec(1), -force_vec(1), theta_p, theta_p, theta_pc, theta_pc, hinge_props.(['c_hinge_' ele_side]), hinge_props.(['c_hinge_' ele_side]), end_disp, end_disp); % Keep residual strength forever
                        
                        % uniaxialMaterial IMKPeakOriented $Mat_Tag $Ke $Up_pos $Upc_pos $Uu_pos $Fy_pos $FmaxFy_pos $FresFy_pos $Up_neg $Upc_neg $Uu_neg $Fy_neg $FmaxFy_neg $FresFy_neg $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $D_pos $D_neg
                        fprintf(fileID,'uniaxialMaterial IMKPeakOriented %i %f %f %f %f %f %f %f %f %f %f %f %f %f 100.0 100.0 100.0 100.0 1.0 1.0 1.0 1.0 1.0 1.0 \n', ele_hinge_id, K0, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength, theta_p, theta_pc, end_disp, f_yield, f_ult_ratio, residual_strength);
               
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

