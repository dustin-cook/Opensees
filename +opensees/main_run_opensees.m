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
if contains(cmdout,'Analysis Failure: Collapse')
    fprintf('Model Reached Collapse Limit \n')
elseif contains(cmdout,'Analysis Failure: Convergence')
    fprintf('Model Experienced a Convergence Failure')
elseif contains(cmdout,'Analysis Failure: Singularity')
    fprintf('Model Experienced a Singularity Failure')
elseif status ~= 0 %(shouldnt get here)
    fprintf('UNEXPECTED OPENSEES FAILURE \n')
else
    fprintf('Model Ran Successfully \n')
end

end

