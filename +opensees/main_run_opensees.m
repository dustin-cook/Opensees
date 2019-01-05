function [ ] = main_run_opensees( opensees_dir )
% Function to trigger the command line to run opensees

%% Run Opensees
command = ['opensees ' opensees_dir filesep 'run_analysis.tcl'];
[status,cmdout] = system(command,'-echo');

% test for analysis failure and terminate Matlab
if contains(cmdout,'Analysis Failure: Convergence') || contains(cmdout,'Analysis Failure: Singularity')
    error('Opensees Failed')
end

end

