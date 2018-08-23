function [ ] = main_run_opensees( output_dir )
% Function to trigger the command line to run opensees

%% Run Opensees
command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
[status,cmdout] = system(command,'-echo');

% test for analysis failure and terminate Matlab
if contains(cmdout,'Analysis Failure:')
    error('Opensees Failed')
end

end

