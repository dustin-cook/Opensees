function [ ] = main_run_opensees( output_dir )
% Function to trigger the command line to run opensees

%% Run Opensees
command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
system(command);

end

