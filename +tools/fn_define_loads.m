function [ ] = fn_define_loads( output_dir, analysis, damp_ratio, node, dt )
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
fprintf(fileID,'constraints Plain \n');
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
for i = 1:length(node.id)
    fprintf(fileID,'  load %d %f 0.0 0.0 0.0 0.0 0.0 \n', node.id(i), node.force(i));
end
fprintf(fileID,'} \n');

%% For Dynamic Analysis
if analysis.type == 3 || analysis.type == 4
    % Define Seismic Excitation Load
    % timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>
    fprintf(fileID,'timeSeries Path 1 -dt %f -filePath %s/%s -factor 386. \n',dt, analysis.eq_dir{1}, analysis.eq_name{1});
    % pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $vel0> <-fact $cFactor>
    fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 \n'); 

    % Define Damping based on first eigen mode
    fprintf(fileID,'set lambda [expr [eigen -fullGenLapack 1]] \n');
    fprintf(fileID,'set omega [expr sqrt($lambda)] \n');
    fprintf(fileID,'set period [expr 6.283185/$omega] \n');
    fprintf(fileID,'puts $period \n');
    fprintf(fileID,'rayleigh 0 0 0 [expr 2*%d/$omega] \n', damp_ratio);  %[expr %d*2*$omega]
    
end

%% Close File
fclose(fileID);

end

