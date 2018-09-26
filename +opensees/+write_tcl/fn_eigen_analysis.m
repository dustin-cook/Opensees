function [ ] = fn_eigen_analysis( output_dir, prim_story_nodes, num_stories, analysis )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Write Recorder File
file_name = [output_dir filesep 'eigen.tcl'];
fileID = fopen(file_name,'w');

% Initail Setup
% Eigen analysis - 10.4 from "Dynamics of Structures" book by Anil Chopra - using equalDOF and very high Ib
fprintf(fileID,'wipeAnalysis \n');		

% Record eigenvectors
fprintf(fileID,'recorder Node -file %s/mode_shape_1.txt -dT %f -node %s -dof 1 3 "eigen 1" \n',output_dir, 1, num2str(prim_story_nodes));
fprintf(fileID,'recorder Node -file %s/mode_shape_2.txt -dT %f -node %s -dof 1 3 "eigen 2" \n',output_dir, 1, num2str(prim_story_nodes));

% Perform Eigen Analysis
if strcmp(analysis.damping,'simple')
    fprintf(fileID,'set numModes %i \n',1);
else
    fprintf(fileID,'set numModes %i \n',min([6,num_stories]));
end
fprintf(fileID,'set lambda [eigen -fullGenLapack $numModes] \n');
fprintf(fileID,'set T {}\n');
fprintf(fileID,'set pi [expr 2.0*asin(1.0)] \n');
fprintf(fileID,'foreach lam $lambda {\n');
fprintf(fileID,'	lappend T [expr (2.0*$pi)/sqrt($lam)]\n');
fprintf(fileID,'}\n');

% Write Output File for periods
fprintf(fileID,'set period_file %s/period.txt \n',output_dir);
fprintf(fileID,'set Periods [open $period_file "w"] \n');
fprintf(fileID,'foreach t $T { \n');
fprintf(fileID,'	puts $Periods " $t" \n');
fprintf(fileID,'} \n');
fprintf(fileID,'close $Periods \n');

% Run Simple one step gravity load to record eigenvectors
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

