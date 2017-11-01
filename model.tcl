#  Define the model (2 dimensions, 3 dof 
model basic -ndm 2 -ndf 3 
# define nodes (inches) 
node 1 0.000000 0.000000 
node 2 0.000000 1000.000000 
# set boundary conditions at each node (3dof) (fix = 1, free = 0) 
fix 1 1 1 1 
fix 2 0 0 0 
# define nodal masses (horizontal) (units?) 
mass 1 0.000000 0. 0. 
mass 2 1.000000 0. 0. 
# Linear Transformation 
geomTransf Linear 1 
# Define Elements (columns and beam) 
# element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag> 
element elasticBeamColumn 1 1 2 9999999999999.000000 60000.000000 5555.555556 1 
