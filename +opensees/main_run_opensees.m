function [ ] = main_run_opensees( opensees_dir, analysis )
% Function to trigger the command line to run opensees

%% Run Opensees
if analysis.summit
    if analysis.opensees_SP
        command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' opensees_dir filesep 'run_analysis.tcl'];
    else
        command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSees ' opensees_dir filesep 'run_analysis.tcl'];
    end
else
    if analysis.opensees_SP
        command = ['openseesSP ' opensees_dir filesep 'run_analysis.tcl'];
    else
        command = ['opensees ' opensees_dir filesep 'run_analysis.tcl'];
    end
end
if analysis.suppress_outputs
    [status,cmdout] = system(command);
else
    [status,cmdout] = system(command,'-echo');
end

% test for analysis failure and terminate Matlab
if contains(cmdout,'Analysis Failure: Convergence') || contains(cmdout,'Analysis Failure: Singularity')
    error('Opensees Failed')
end

end

