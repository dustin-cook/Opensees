function [ ] = vezna_sdof( period, analysis, dt, num_steps )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% define I for this run
I = (5.18 * 432^3) / (3 * period^2 * 3225 );

%% Write TCL file
file_name = 'vezna_sdof.tcl';
fileID = fopen(file_name,'w');

fprintf(fileID,'wipe; \n');
fprintf(fileID,'model basic -ndm 2 -ndf 3; \n');
fprintf(fileID,'file mkdir data; \n');
fprintf(fileID,'node 1 0. 0.; \n');
fprintf(fileID,'node 2 0. 432. \n');
fprintf(fileID,'fix 1 1 1 1; \n');
fprintf(fileID,'mass 2 5.18 0. 0.; \n');
fprintf(fileID,'geomTransf Linear 1; \n');
fprintf(fileID,'element elasticBeamColumn 1 1 2 3600 3225 %f 1; \n',I);
% fprintf(fileID,'recorder Node -file Data/DFree.out -time -node 2 -dof 1 2 3 disp; \n');
% fprintf(fileID,'recorder Node -file Data/RBase.out -time -node 1 -dof 1 2 3 reaction; \n');
% fprintf(fileID,'recorder Drift -file Data/Drift.out -time -iNode 1 -jNode 2 -dof 1  -perpDirn 2 ; \n');
% fprintf(fileID,'recorder Element -file Data/FCol.out -time -ele 1 force; \n');
fprintf(fileID,'recorder Node -file vezna_accel.txt -time -node 2 -dof 1 accel; \n');
fprintf(fileID,'timeSeries Linear 1 \n');
fprintf(fileID,'pattern Plain 1 1 { \n');
fprintf(fileID,'   load 2 0. -2000. 0.; \n');
fprintf(fileID,'} \n');
fprintf(fileID,'constraints Plain; \n');
fprintf(fileID,'numberer Plain; \n');
fprintf(fileID,'system BandGeneral; \n');
fprintf(fileID,'algorithm Linear; \n');
fprintf(fileID,'integrator LoadControl 0.1; \n');
fprintf(fileID,'analysis Static \n');
fprintf(fileID,'analyze 10; \n');
fprintf(fileID,'loadConst -time 0.0; \n');
fprintf(fileID,'set G 386 \n');
% fprintf(fileID,'timeSeries Path 2 -dt %f -filePath %s/%s -factor $G; \n',dt, analysis.eq_dir, analysis.eq_name);
fprintf(fileID,'timeSeries Path 2 -dt %f -filePath A10000.tcl -factor $G; \n',dt);
fprintf(fileID,'pattern UniformExcitation 2 1 -accel 2; \n');
fprintf(fileID,'set freq [expr [eigen -fullGenLapack 1]**0.5] \n');
fprintf(fileID,'set dampRatio 0.02 \n');
fprintf(fileID,'rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq] \n');
% fprintf(fileID,'recorder display "Displaced shape" 10 10 500 500 -wipe \n');
% fprintf(fileID,'prp 200. 50. 1; \n');
% fprintf(fileID,'vup  0  1 0; \n');
% fprintf(fileID,'vpn  0  0 1; \n');
% fprintf(fileID,'display 1 5 40  \n');
fprintf(fileID,'wipeAnalysis; \n');
fprintf(fileID,'constraints Plain; \n');
fprintf(fileID,'numberer Plain; \n');
fprintf(fileID,'system BandGeneral; \n');
fprintf(fileID,'algorithm Linear \n');
fprintf(fileID,'integrator Newmark 0.5 0.25 ; \n');
fprintf(fileID,'analysis Transient; \n');
fprintf(fileID,'analyze %d %f; \n',num_steps,analysis.time_step);
fprintf(fileID,'puts "Done!" \n');
fprintf(fileID,'wipe \n');

% Close File
fclose(fileID);

end

