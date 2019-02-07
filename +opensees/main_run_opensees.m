function [ ] = main_run_opensees( opensees_dir, analysis )
% Function to trigger the command line to run opensees

%% Run Opensees
if analysis.summit
    command = ['mpirun -np $SLURM_NTASKS /projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' opensees_dir filesep 'run_analysis.tcl'];
elseif analysis.opensees_SP
    command = ['openseesSP ' opensees_dir filesep 'run_analysis.tcl'];
else
    command = ['opensees ' opensees_dir filesep 'run_analysis.tcl'];
end
[status,cmdout] = system(command,'-echo');

% test for analysis failure and terminate Matlab
if contains(cmdout,'Analysis Failure: Convergence') || contains(cmdout,'Analysis Failure: Singularity')
    error('Opensees Failed')
end

end

