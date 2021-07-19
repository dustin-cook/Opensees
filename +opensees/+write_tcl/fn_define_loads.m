function [ ] = fn_define_loads( write_dir, analysis, node, dimension, story, element, joint, ground_motion, model )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup 
% Import Packages
import asce_7.*

% Define site hazard design values (pass into function instead)
site.S1 = 0.9/1.4; %fv from curts site assuming site class C
site.Sds = 1;
site.Sd1 = 0.6;

% Write Loads File
file_name = [write_dir filesep 'loads.tcl'];
fileID = fopen(file_name,'w');

%% Define Gravity Loads (node id, axial, shear, moment)
fprintf(fileID,'puts "Defining Loads ..." \n');
fprintf(fileID,'pattern Plain 1 Linear {  \n');
for i = 1:height(element) 
    if strcmp(element.type{i},'beam')
        element.w(i,1) = element.gravity_load(i)/element.length(i);
        % eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wy <$Wx>
        % eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wy $Wz <$Wx>
        fprintf(fileID,'   eleLoad -ele %i -type -beamUniform -%d 0.0 \n', element.id(i), element.w(i));
    else % walls and columns
        if strcmp(dimension,'2D')
            ele_grav_load_1 = element.gravity_load_1(i);
            ele_grav_load_2 = element.gravity_load_2(i);
            % eleLoad -range $eleTag1 $eleTag2 -type -beamPoint $Py $xL <$Px>
            fprintf(fileID,'   eleLoad -ele %d -type -beamPoint 0.0 0.0 -%d \n', element.id(i), ele_grav_load_1); 
            fprintf(fileID,'   eleLoad -ele %d -type -beamPoint 0.0 1.0 -%d \n', element.id(i), ele_grav_load_2);
        elseif strcmp(dimension,'3D')
            % eleLoad -range $eleTag1 $eleTag2 -type -beamPoint $Py $Pz $xL <$Px>
            fprintf(fileID,'   eleLoad -ele %d -type -beamPoint 0.0 0.0 0.0 -%d \n', element.id(i), ele_grav_load_1); 
            fprintf(fileID,'   eleLoad -ele %d -type -beamPoint 0.0 0.0 1.0 -%d \n', element.id(i), ele_grav_load_2); 
        end
    end
end

% Add distributed loads to joint section of beams (just to archetype
% models)
if analysis.model_type == 3
    for i = 1:height(joint)
    % beam left
    if joint.beam_left(i) > 0
        joint_ele_id = joint.id(i)*1000000+1;
        dist_load = element.w(element.id == joint.beam_left(i));
        fprintf(fileID,'   eleLoad -ele %i -type -beamUniform -%d 0.0 \n', joint_ele_id, dist_load);
    end
    % beam right
    if joint.beam_right(i) > 0
        joint_ele_id = joint.id(i)*1000000+2;
        dist_load = element.w(element.id == joint.beam_right(i));
        fprintf(fileID,'   eleLoad -ele %i -type -beamUniform -%d 0.0 \n', joint_ele_id, dist_load);
    end
    
    end
end


% for i = 1:height(node) 
%     % eleLoad -ele $eleTag1 <$eleTag2 ....> -type -beamUniform $Wy $Wz <$Wx>
%     if node.mass(i) > 0
%         fprintf(fileID,'   load %i 0.0 %f 0.0 \n', node.id(i), -386*node.mass(i));
%     end
% end
fprintf(fileID,'} \n');

%% Write Gravity System Analysis
fprintf(fileID,'constraints Transformation \n');
fprintf(fileID,'numberer RCM \n'); % renumber dof's to minimize band-width (optimization)
if analysis.opensees_SP
    fprintf(fileID,'system Mumps ICNTL 100 \n'); % Use Mumps for OpenseesSP
else
    fprintf(fileID,'system BandGeneral \n'); % how to store and solve the system of equations in the analysis
end
fprintf(fileID,'test NormDispIncr 1.0e-5 1000 \n');
% fprintf(fileID,'test EnergyIncr 1.0e-6 50 \n'); % determine if convergence has been achieved at the end of an iteration step
% fprintf(fileID,'algorithm Linear \n');
fprintf(fileID,'algorithm KrylovNewton \n');
fprintf(fileID,'integrator LoadControl 0.1 \n');
fprintf(fileID,'analysis Static \n');
fprintf(fileID,'set ok [analyze 10] \n');
fprintf(fileID,'if {$ok != 0} { \n');
fprintf(fileID,'puts "Analysis Failure: Gravity Load Failure" \n');
fprintf(fileID,'wipe \n');
fprintf(fileID,'exit \n');
fprintf(fileID,'} \n');
if ~analysis.suppress_outputs
    fprintf(fileID,'puts "Gravity Load Complete" \n');
end
fprintf(fileID,'loadConst -time 0.0 \n');

%% Static Lateral Loading
if analysis.type == 2 || analysis.type == 3 || analysis.type == 4
    % equivalent lateral force vertical distribution
    w = story.seismic_wt;
    h = story.story_ht + story.y_start;
    if ~isfield(analysis,'run_drifts')
        analysis.run_drifts = 0; %default not to run if it hasn't been defined
    end
    [ ~, ~, story.lateral_force, ~ ] = fn_equivalent_lateral_force( w, h, model.T1_x, site.S1, site.Sds, site.Sd1, model.ie, analysis.run_drifts );
    
    % define nodes to push
    force_nodes = [];
    node.lateral_load = zeros(height(node),1);
    for s = 1:height(story)
        if story.id(s) > 0
            force_nodes_this_story = node.id(node.story == story.id(s) & node.mass > 0);
            node.lateral_load(ismember(node.id,force_nodes_this_story)) = analysis.eq_lat_load_factor * story.lateral_force(s)/length(force_nodes_this_story);
            force_nodes = [force_nodes;force_nodes_this_story];
        end
    end
    
    % Define Static Lateral Load Pattern
    fprintf(fileID,'pattern Plain 2 Linear { \n');
    for i = 1:length(force_nodes)
        node_id = force_nodes(i);
        if strcmp(dimension,'2D')
            if strcmp(analysis.pushover_direction,'x')
                fprintf(fileID,'  load %d %f 0.0 0.0 \n', node_id, node.lateral_load(node.id == node_id));
            elseif strcmp(analysis.pushover_direction,'-x')
                fprintf(fileID,'  load %d %f 0.0 0.0 \n', node_id, -node.lateral_load(node.id == node_id));
            end
        elseif strcmp(dimension,'3D')
            if strcmp(analysis.pushover_direction,'x')
                fprintf(fileID,'  load %d %f 0.0 0.0 0.0 0.0 0.0 \n', node_id, node.lateral_load(node.id == node_id));
            elseif strcmp(analysis.pushover_direction,'-x')
                fprintf(fileID,'  load %d %f 0.0 0.0 0.0 0.0 0.0 \n', node_id, -node.lateral_load(node.id == node_id));
            elseif strcmp(analysis.pushover_direction,'z')
                fprintf(fileID,'  load %d 0.0 0.0 %f 0.0 0.0 0.0 \n', node_id, node.lateral_load(node.id == node_id));
            elseif strcmp(analysis.pushover_direction,'-z')
                fprintf(fileID,'  load %d 0.0 0.0 %f 0.0 0.0 0.0 \n', node_id, -node.lateral_load(node.id == node_id));
            end
        end
    end
    fprintf(fileID,'} \n');
end

%% Dynamic Analysis
if analysis.type == 1
    scale_factor = 386;%*analysis.ground_motion_scale_factor; % g's to in per s times scale factor
    % Define Seismic Excitation Load
    % timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>
    % pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $vel0> <-fact $cFactor>
    if isfield(ground_motion,'x')
        fprintf(fileID,'set dt %f \n',ground_motion.x.eq_dt);
%         fprintf(fileID,'puts "EQ X dt = $dt" \n');
        fprintf(fileID,'timeSeries Path 1 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.x.eq_dir{1}, ground_motion.x.eq_name{1}, scale_factor);
        fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 -fact %f \n',analysis.ground_motion_scale_factor); 
    end
    if isfield(ground_motion,'z') && strcmp(dimension,'3D')
        fprintf(fileID,'set dt %f \n',ground_motion.z.eq_dt);
%         fprintf(fileID,'puts "EQ Z dt = $dt" \n');
        fprintf(fileID,'timeSeries Path 2 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.z.eq_dir{1}, ground_motion.z.eq_name{1}, scale_factor);
        fprintf(fileID,'pattern UniformExcitation 4 3 -accel 2 -fact %f \n',analysis.ground_motion_scale_factor); 
    end
    if isfield(ground_motion,'y')
        fprintf(fileID,'set dt %f \n',ground_motion.y.eq_dt);
%         fprintf(fileID,'puts "EQ Y dt = $dt" \n');
        fprintf(fileID,'timeSeries Path 3 -dt $dt -filePath %s/%s -factor %f \n', ground_motion.y.eq_dir{1}, ground_motion.y.eq_name{1}, scale_factor);
        fprintf(fileID,'pattern UniformExcitation 5 2 -accel 3 -fact %f \n',analysis.ground_motion_scale_factor); 
    end
end


% Define Damping based on eigen modes
if analysis.type == 1
    if ~analysis.suppress_outputs
        fprintf(fileID,'puts "Running Eigen and Defining Damping" \n');
    end
    if  strcmp(analysis.damping,'simple')
        fprintf(fileID,'set lambda [eigen 1] \n');
        fprintf(fileID,'set pi [expr 2.0*asin(1.0)] \n');
        fprintf(fileID,'set omega [expr sqrt($lambda)] \n');
        fprintf(fileID,'set period [expr 2*$pi/sqrt($lambda)] \n');
        fprintf(fileID,'puts $period \n');
        fprintf(fileID,'set alpha [expr %d*$omega] \n', analysis.damp_ratio);
        fprintf(fileID,'set beta [expr %d/$omega] \n', analysis.damp_ratio);
        fprintf(fileID,'rayleigh $alpha 0.0 $beta 0.0 \n'); 
    else
        fprintf(fileID,'set lambda [eigen %i] \n', 3);
        fprintf(fileID,'set pi [expr 2.0*asin(1.0)] \n');
        fprintf(fileID,'set i 0 \n');
        fprintf(fileID,'foreach lam $lambda {\n');
        fprintf(fileID,'    set i [expr $i+1] \n');
        fprintf(fileID,'    puts $lam \n');
        fprintf(fileID,'	set omega($i) [expr sqrt($lam)]\n');
        fprintf(fileID,'	set period($i) [expr 2*$pi/sqrt($lam)]\n');
        fprintf(fileID,'}\n');
        fprintf(fileID,'puts $period(1) \n');
        fprintf(fileID,'puts $period(2) \n');
        fprintf(fileID,'puts $period(3) \n');
        if  strcmp(analysis.damping,'modal')
            fprintf(fileID,'set zeta %f\n',0.002);		% percentage of critical damping
        else
            fprintf(fileID,'set zeta %f\n',analysis.damp_ratio);		% percentage of critical damping
        end
        fprintf(fileID,'set a0 [expr $zeta*2.0*$omega(1)*$omega(2)/($omega(1) + $omega(2))]\n');	% mass damping coefficient based on first and third modes
        if analysis.nonlinear == 0
            fprintf(fileID,'set a1 [expr $zeta*2.0/($omega(1) + $omega(2))]\n'); % stiffness damping coefficient based on first and third modes
        else
            % Modify Stiffness Proportional Coefficient for Nonlinear hinge model according to Ibbara 2005
            stiffness_mod = (analysis.hinge_stiff_mod+1)/analysis.hinge_stiff_mod;
            fprintf(fileID,'set a1 [expr %f*$zeta*2.0/($omega(1) + $omega(2))]\n',stiffness_mod); % stiffness damping coefficient based on first and third modes
        end

        % Set Damping
        if strcmp(analysis.damping,'rayleigh')
            % region $regTag <-ele ($ele1 $ele2 ...)> <-eleOnly ($ele1 $ele2 ...)> <-eleRange $startEle $endEle> <-eleOnlyRange $startEle $endEle> <-node ($node1 $node2 ...)> <-nodeOnly ($node1 $node2 ...)> <-nodeRange $startNode $endNode> <-nodeOnlyRange $startNode $endNode> <-node all> <-rayleigh $alphaM $betaK $betaKinit $betaKcomm>
    %         rayleigh $alphaM $betaK $betaKinit $betaKcomm
    %         fprintf(fileID,'rayleigh $a0 $a1 0.0 0.0 \n');
    %         Region Commands do not work (ie they apply no damping)
            if analysis.nonlinear == 0
                fprintf(fileID,'region 1 -ele %s -rayleigh 0.0 $a1 0.0 0.0 \n', num2str(element.id')); % Assign Stiffnes Proportional Damping to the elastic elements
            else
                fprintf(fileID,'region 1 -ele %s -rayleigh 0.0 $a1 0.0 0.0 \n', num2str(element.id(element.rigid == 0)')); % Assign Stiffnes Proportional Damping to the elastic elements
            end
            fprintf(fileID,'region 2 -node %s -rayleigh $a0 0.0 0.0 0.0 \n', num2str(node.id(node.mass > 0)')); % Assign Mass Proportional Damping to the whole model (only triggers where there is mass)
        elseif strcmp(analysis.damping,'modal')
                fprintf(fileID,'modalDamping %f \n',analysis.damp_ratio);
                fprintf(fileID,'rayleigh 0.0 $a1 0.0 0.0 \n'); % Just a little bit of stiffness proportional everywhere
        else
            error('Damping Type Not Recognized')
        end
    end
end

if ~analysis.suppress_outputs
    fprintf(fileID,'puts "Define Load Complete" \n');
end

%% Close File
fclose(fileID);

end

