function [ ] = fn_eigen_analysis( write_dir, prim_story_nodes, num_stories, analysis, dims )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define File Type to Write to
if analysis.write_xml
    file_type = '-xml';
    file_ext = 'xml';
else
    file_type = '-file';
    file_ext = 'txt';
end

% Write Eigen File
file_name = [write_dir filesep 'eigen.tcl'];
fileID = fopen(file_name,'w');

fprintf(fileID,'puts "Running Eigen ..." \n');

% Initail Setup
if strcmp(analysis.damping,'simple')
    num_modes = 2;
else
    if strcmp(dims,'3D')
        num_modes = min([6,num_stories+1]);
    else
        num_modes = min([6,num_stories]);
    end
end

% Eigen analysis - 10.4 from "Dynamics of Structures" book by Anil Chopra - using equalDOF and very high Ib
fprintf(fileID,'wipeAnalysis \n');		

% Record eigenvectors
for i = 1:num_modes
    fprintf(fileID,'recorder Node %s %s/mode_shape_%i.%s -dT %f -node %s -dof 1 3 "eigen %s" \n', file_type, write_dir, i, file_ext, 1, num2str(prim_story_nodes), num2str(i));
end

% Perform Eigen Analysis
fprintf(fileID,'set numModes %i \n',num_modes);
if strcmp(analysis.damping,'simple')
    fprintf(fileID,'set lambda [eigen -fullGenLapack $numModes ] \n');
else
    fprintf(fileID,'set lambda [eigen $numModes ] \n');
end
fprintf(fileID,'set T {}\n');
fprintf(fileID,'set pi [expr 2.0*asin(1.0)] \n');
fprintf(fileID,'foreach lam $lambda {\n');
fprintf(fileID,'	lappend T [expr (2.0*$pi)/sqrt($lam)]\n');
fprintf(fileID,'}\n');

% Write Output File for periods
fprintf(fileID,'set period_file %s/period.txt \n',write_dir);
fprintf(fileID,'set Periods [open $period_file "w"] \n');
fprintf(fileID,'foreach t $T { \n');
fprintf(fileID,'	puts $Periods " $t" \n');
fprintf(fileID,'} \n');
fprintf(fileID,'close $Periods \n');

% Run Simple one step gravity load to record eigenvectors
fprintf(fileID,'integrator LoadControl 0 1 0 0 \n');
fprintf(fileID,'test EnergyIncr 1.0e-10 100 0 \n');
fprintf(fileID,'algorithm KrylovNewton \n');
fprintf(fileID,'numberer RCM \n');
fprintf(fileID,'constraints Transformation \n');
if analysis.opensees_SP
    fprintf(fileID,'system Mumps ICNTL 100 \n'); % Use Mumps for OpenseesSP
else
    fprintf(fileID,'system ProfileSPD \n');
end
fprintf(fileID,'analysis Static \n');
fprintf(fileID,'set res [analyze 1] \n');
fprintf(fileID,'if {$res < 0} { \n');
fprintf(fileID,'    puts "Modal analysis failed" \n');
fprintf(fileID,'} \n');

% Remove recorder
fprintf(fileID,'remove recorders \n');

if ~analysis.suppress_outputs
    fprintf(fileID,'puts "Eigen Complete" \n');
end

% Close File
fclose(fileID);
end

