function [ ] = fn_define_loads( output_dir, analysis, damp_ratio, node, ground_motion )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Write Loads File
file_name = [output_dir filesep 'loads.tcl'];
fileID = fopen(file_name,'w');

%% Define Gravity Loads (node id, axial, shear, moment)
fprintf(fileID,'pattern Plain 1 Linear {  \n');
for i = 1:length(node.id)
    fprintf(fileID,'   load %d 0.0 -%f 0.0 0.0 0.0 0.0 \n', node.id(i), node.weight(i));
end
fprintf(fileID,'} \n');

% Write Gravity System Analysis
fprintf(fileID,'constraints Transformation \n');
fprintf(fileID,'numberer RCM \n'); % renumber dof's to minimize band-width (optimization)
fprintf(fileID,'system BandGeneral \n'); % how to store and solve the system of equations in the analysis
fprintf(fileID,'test EnergyIncr 0.00000001 6 \n'); % determine if convergence has been achieved at the end of an iteration step
fprintf(fileID,'algorithm Newton \n');
fprintf(fileID,'integrator LoadControl 0.1 \n');
fprintf(fileID,'analysis Static	 \n');
fprintf(fileID,'analyze 10 \n');
fprintf(fileID,'loadConst -time 0.0 \n');

%% Define Static Lateral Load Patter
fprintf(fileID,'pattern Plain 2 Linear { \n');
node_force = 0; % Just turn off for now
for i = 1:length(node.id)
    fprintf(fileID,'  load %d %f 0.0 0.0 0.0 0.0 0.0 \n', node.id(i), node_force);
end
fprintf(fileID,'} \n');

%% For Dynamic Analysis
if analysis.type == 3 || analysis.type == 4
    % Define Seismic Excitation Load
    % timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>
    % pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $vel0> <-fact $cFactor>
    if ~strcmp(ground_motion.eq_name_x,'')
        fprintf(fileID,'timeSeries Path 1 -dt %f -filePath %s/%s -factor 386. \n',ground_motion.eq_dt, ground_motion.eq_dir{1}, ground_motion.eq_name_x{1});
        fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 -fact %f \n',ground_motion.x_ratio); 
    end
    if ~strcmp(ground_motion.eq_name_z,'')
        fprintf(fileID,'timeSeries Path 2 -dt %f -filePath %s/%s -factor 386. \n',ground_motion.eq_dt, ground_motion.eq_dir{1}, ground_motion.eq_name_z{1});
        fprintf(fileID,'pattern UniformExcitation 4 3 -accel 2 -fact %f \n',ground_motion.z_ratio); 
    end
    if ~strcmp(ground_motion.eq_name_y,'')
        fprintf(fileID,'timeSeries Path 3 -dt %f -filePath %s/%s -factor 386. \n',ground_motion.eq_dt, ground_motion.eq_dir{1}, ground_motion.eq_name_y{1});
        fprintf(fileID,'pattern UniformExcitation 5 2 -accel 3 -fact %f \n',ground_motion.y_ratio); 
    end

    % Define Damping based on first eigen mode
    fprintf(fileID,'set lambda [eigen -fullGenLapack 3] \n');
    fprintf(fileID,'set pi 3.141593\n');
    fprintf(fileID,'set i 0 \n');
    fprintf(fileID,'foreach lam $lambda {\n');
    fprintf(fileID,'    set i [expr $i+1] \n');
    fprintf(fileID,'	set omega($i) [expr sqrt($lam)]\n');
    fprintf(fileID,'	set period($i) [expr 2*$pi/sqrt($lam)]\n');
    fprintf(fileID,'}\n');
    fprintf(fileID,'puts $period(1) \n');
    fprintf(fileID,'puts $period(3) \n');
    fprintf(fileID,'set alpha [expr 2*%d*(1-$omega(1))/(1/$omega(1) - $omega(1)/($omega(3)*$omega(3)))]\n', damp_ratio);
    fprintf(fileID,'set beta [expr 2*%d - $alpha/($omega(3)*$omega(3))]\n', damp_ratio);
    fprintf(fileID,'rayleigh $alpha 0 $beta 0 \n');  
end

%% Close File
fclose(fileID);

end

