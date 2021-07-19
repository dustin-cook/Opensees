function [ ] = fn_define_eignen_run_script( write_dir )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import opensees.write_tcl.*

%% Write Loads File
file_name = [write_dir filesep 'run_analysis.tcl'];
fileID = fopen(file_name,'w');

%% Initial Analysis Setup
fprintf(fileID,'wipe \n');
fprintf(fileID,'source %s/model.tcl \n', write_dir);
fprintf(fileID,'source %s/eigen.tcl \n', write_dir);

% Close File
fclose(fileID);

end

