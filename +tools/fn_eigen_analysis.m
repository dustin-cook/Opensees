function [ ] = fn_eigen_analysis( output_dir, time_step, nodes )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Write Recorder File
file_name = [output_dir filesep 'eigen.tcl'];
fileID = fopen(file_name,'w');

% Initail Setup
fprintf(fileID,'# Eigen analysis - 10.4 from "Dynamics of Structures" book by Anil Chopra - using equalDOF and very high Ib \n');
fprintf(fileID,'wipeAnalysis \n');		

% Record eigenvectors
fprintf(fileID,'# record eigenvectors \n');
fprintf(fileID,'recorder Node -file A%s/mode_shape_1.txt -dT %f -node %s -dof 1 "eigen 1" \n',output_dir, 2*time_step, num2str(nodes'));

% Perform Eigen Analysis
fprintf(fileID,'# perform eigen analysis\n');
fprintf(fileID,'set numModes %d \n',6);
fprintf(fileID,'set lambda [eigen -fullGenLapack $numModes] \n');
fprintf(fileID,'set T {}\n');
fprintf(fileID,'set pi 3.141593\n');
fprintf(fileID,'foreach lam $lambda {\n');
fprintf(fileID,'	lappend T [expr (2*$pi)/sqrt($lam)]\n');
fprintf(fileID,'}\n');

% Write Output File for periods
fprintf(fileID,'# write the output file cosisting of periods \n');
fprintf(fileID,'set period_file %s/period.txt \n',output_dir);
fprintf(fileID,'set Periods [open $period_file "w"] \n');
fprintf(fileID,'foreach t $T { \n');
fprintf(fileID,'	puts $Periods " $t" \n');
fprintf(fileID,'} \n');
fprintf(fileID,'close $Periods \n');

% Run Simple one step gravity load to record eigenvectors
fprintf(fileID,'# Run a one step gravity load with no loading (to record eigenvectors)\n');
fprintf(fileID,'integrator LoadControl 0 1 0 0 \n');
fprintf(fileID,'test EnergyIncr 1.0e-10 100 0 \n');
fprintf(fileID,'algorithm Newton \n');
fprintf(fileID,'numberer RCM \n');
fprintf(fileID,'constraints Transformation \n');
fprintf(fileID,'system ProfileSPD \n');
fprintf(fileID,'analysis Static \n');
fprintf(fileID,'set res [analyze 1] \n');
fprintf(fileID,'if {$res < 0} { \n');
fprintf(fileID,'    puts "Modal analysis failed" \n');
fprintf(fileID,'} \n');

% Remove recorder
fprintf(fileID,'remove recorders \n');
end

