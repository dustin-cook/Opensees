% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup
import tools.*

%% Load Analysis and Model parameters
analysis = readtable(['inputs' filesep 'analysis.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
% Create Outputs Directory
output_dir = ['outputs/' model.name{1}];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Define EQ ground motion
eqs = dir([analysis.eq_dir{1} filesep '*eq_*']);
num_eqs = length(eqs);

% Main Opensees Analysis
for i = 1:1%length(eqs)
    % EQ this run
    eq = load([analysis.eq_dir{1} filesep analysis.eq_name{1}]);
        
    %% Create Model Databases
    [ node, element, story_ht, num_node_story ] = fn_model_table( model );

    %% Write TCL file
    fn_build_model( output_dir, node, element )
    fn_define_recorders( output_dir, analysis.type, node.id, element.id )
    fn_define_loads( output_dir, analysis, model.damp_ratio, node, analysis.eq_dt )
    fn_define_analysis( output_dir, analysis, node.id, eq, analysis.eq_dt )
    fn_eigen_analysis( output_dir, analysis.time_step, node.id )

    %% Run Opensees
    command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
    system(command);
    
    %% Save workspace data
    save([output_dir filesep 'analysis_data'])
    
end
toc

