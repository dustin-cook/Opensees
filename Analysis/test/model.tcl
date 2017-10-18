#  Define the model (2 dimensions, 3 dof 
model basic -ndm 2 -ndf 3 
# define nodes (inches) 
node 1 0 0 
node 2 240 0 
node 3 0 120 
node 4 240 120 
# set boundary conditions at each node (3dof) (fix = 1, free = 0) 
fix 1 1 1 0 
fix 2 1 1 0 
fix 3 0 0 0 
fix 4 0 0 0 
# define nodal masses (horizontal) (units?) 
mass 3 0 0. 0. 
mass 4 0 0. 0. 
# Linear Transformation 
geomTransf Linear 1 
# Define Elements (columns and beam) 
# element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag> 
element elasticBeamColumn 1 1 3 9999999999 9999999999 9999999999 1 
element elasticBeamColumn 2 2 4 9999999999 9999999999 9999999999 1 
element elasticBeamColumn 3 3 4 9999999999 9999999999 9999999999 1 
