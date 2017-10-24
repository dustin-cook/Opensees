function [ ] = fn_define_loads( analysis, damp_ratio, node, dt )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Write Loads File
file_name = ['Analysis' filesep analysis.id filesep 'loads.tcl'];
fileID = fopen(file_name,'w');

% Define Gravity Loads (node id, axial, shear, moment)
fprintf(fileID,'# Define Gravity Loads (node id, axial, shear, moment) \n');
fprintf(fileID,'pattern Plain 1 Linear {  \n');
for i = 1:length(node.id)
    fprintf(fileID,'   load %d 0.0 %d 0.0 \n', node.id(i), node.weight(i));
end
fprintf(fileID,'} \n');

% Define Static Lateral Load Patter
fprintf(fileID,'# Define Load Pattern \n');
fprintf(fileID,'pattern Plain 2 Linear { \n');
for i = 1:length(node.id)
    fprintf(fileID,'  load %d %d 0.0 0.0 \n', node.id(i), node.force(i));
end
fprintf(fileID,'} \n');

% For Dynamic Analysis
if analysis.type == 3 || analysis.type == 4
    % Define Seismic Excitation Load
    fprintf(fileID,'# Define Seismic Excitation Pattern \n');
    fprintf(fileID,'timeSeries Path 1 -dt %d -filePath %s/%s -factor 386 \n',dt, analysis.eq_dir, analysis.eq_name); % timeSeries Path $tag -dt $dt -filePath $filePath <-factor $cFactor> <-useLast> <-prependZero> <-startTime $tStart>
    fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 \n'); % pattern UniformExcitation $patternTag $dir -accel $tsTag <-vel0 $vel0> <-fact $cFactor>

    % Define Damping
    fprintf(fileID,'# set damping based on first eigen mode \n');
    fprintf(fileID,'set freq [expr [eigen -fullGenLapack 1]**0.5] \n');
    fprintf(fileID,'set period [expr 1/$freq] \n');
    fprintf(fileID,'  rayleigh 0. 0. 0. [expr 2*%d*$period] \n', damp_ratio);
    fprintf(fileID,'puts $period \n');
end

% Close File
fclose(fileID);

end

