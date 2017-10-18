function [ ] = fn_build_model( story_ht_in, bay_width_in, foundation_fix, story_mass, E, A, I, analysis_id )
%UNTITLED6 Summary of this function goes here
% ASSUMPTIONS:
% 2d linear model single story, 1 bay, uniform members

%% Write TCL file
file_name = ['Analysis' filesep analysis_id filesep 'model.tcl'];
fileID = fopen(file_name,'w');

% Define the model (2 dimensions, 3 dof
fprintf(fileID,'#  Define the model (2 dimensions, 3 dof \n');
fprintf(fileID,'model basic -ndm 2 -ndf 3 \n');

% define nodes (inches)
fprintf(fileID,'# define nodes (inches) \n');
fprintf(fileID,'node 1 0 0 \n');
fprintf(fileID,'node 2 %d 0 \n',bay_width_in);
fprintf(fileID,'node 3 0 %d \n',story_ht_in);
fprintf(fileID,'node 4 %d %d \n',bay_width_in,story_ht_in);

% set boundary conditions at each node (3dof) (fix = 1, free = 0)
fprintf(fileID,'# set boundary conditions at each node (3dof) (fix = 1, free = 0) \n');
fprintf(fileID,'fix 1 %d %d %d \n',foundation_fix(1),foundation_fix(2),foundation_fix(3));
fprintf(fileID,'fix 2 %d %d %d \n',foundation_fix(1),foundation_fix(2),foundation_fix(3));
fprintf(fileID,'fix 3 0 0 0 \n');
fprintf(fileID,'fix 4 0 0 0 \n');

% define nodal masses (horizontal) (units?)
fprintf(fileID,'# define nodal masses (horizontal) (units?) \n');
fprintf(fileID,'mass 3 %d 0. 0. \n',story_mass/2);
fprintf(fileID,'mass 4 %d 0. 0. \n',story_mass/2);

% Linear Transformation
fprintf(fileID,'# Linear Transformation \n');
fprintf(fileID,'geomTransf Linear 1 \n');

% Define Elements (columns and beam)
% element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag>
fprintf(fileID,'# Define Elements (columns and beam) \n');
fprintf(fileID,'# element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag> \n');
fprintf(fileID,'element elasticBeamColumn 1 1 3 %d %d %d 1 \n',A,E,I);
fprintf(fileID,'element elasticBeamColumn 2 2 4 %d %d %d 1 \n',A,E,I);
fprintf(fileID,'element elasticBeamColumn 3 3 4 %d %d %d 1 \n',A,E,I);

% Close File
fclose(fileID);

end

