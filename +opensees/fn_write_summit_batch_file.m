function [ ] = fn_write_summit_batch_file( write_dir, analysis_name, model_name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Write Eigen File
file_name = [write_dir filesep 'run_summit.sh'];
fileID = fopen(file_name,'w');

% Write Batch Commands
fprintf(fileID,'#!/bin/sh \n');
fprintf(fileID,'#SBATCH --partition=shas # Summit partition \n');	
fprintf(fileID,'#SBATCH --qos=normal # Summit qos \n');	
fprintf(fileID,'#SBATCH --time=24:00:00 # Max wall time \n');	
fprintf(fileID,'#SBATCH --ntasks=1 # Number of tasks per job \n');	
fprintf(fileID,'#SBATCH --nodes=1 # Number of nodes per job \n');	
fprintf(fileID,'#SBATCH --job-name=%s # Job submission name \n', analysis_name);	
fprintf(fileID,'#SBATCH --output=opensees.%s.out # Output file name with Job ID \n', '%j');	
fprintf(fileID,'#SBATCH --mail-type=END # Email user when job finishes \n');	
fprintf(fileID,'#SBATCH --mail-user=dustin.cook@colorado.edu # Email address of user \n');
fprintf(fileID,'\n');

% Documentation
fprintf(fileID,'# Written by:	 Dustin Cook \n');	
fprintf(fileID,'# Date:		     19 December 2018 \n');	
fprintf(fileID,'# Purpose: 	   This script submits an opensees job to the Slurm job scheduler \n');
fprintf(fileID,'\n');

% Load Modules
fprintf(fileID,'# purge all existing modules \n');	
fprintf(fileID,'module purge \n');
fprintf(fileID,'\n');
fprintf(fileID,'# load Opensees and Matlab module \n');	
fprintf(fileID,'module load matlab/2018b \n'); % Load Matlab	
fprintf(fileID,'module load intel/17.4 impi/17.3 hdf5/1.10.1 mkl/17.3 \n');	% Load Opensees
fprintf(fileID,'\n');

% Defined the Job Directory
fprintf(fileID,'# Defined the Job Directory \n');	
fprintf(fileID,'cd /scratch/summit/duco1061/%s/%s \n', model_name, analysis_name);	
fprintf(fileID,'\n');

% Run Command
fprintf(fileID,'# Run Opensees \n');	
fprintf(fileID,'mpirun -np 1 /projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP run_analysis.tcl \n');	

% Close File
fclose(fileID);

end

