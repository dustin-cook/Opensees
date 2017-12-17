#  Define the model (2 dimensions, 3 dof 
model basic -ndm 2 -ndf 3 
# define nodes (inches) 
node 1 0.000000 0.000000 
node 2 0.000000 198.000000 
node 3 0.000000 360.000000 
node 4 0.000000 522.000000 
node 5 0.000000 684.000000 
node 6 0.000000 846.000000 
node 7 0.000000 1004.400000 
# set boundary conditions at each node (3dof) (fix = 1, free = 0) 
fix 1 1 1 1 
fix 2 0 0 0 
fix 3 0 0 0 
fix 4 0 0 0 
fix 5 0 0 0 
fix 6 0 0 0 
fix 7 0 0 0 
# define nodal masses (horizontal) (units?) 
mass 1 0.000000 0. 0. 
mass 2 5.699482 0. 0. 
mass 3 6.023316 0. 0. 
mass 4 6.023316 0. 0. 
mass 5 6.023316 0. 0. 
mass 6 6.023316 0. 0. 
mass 7 4.922280 0. 0. 
# Linear Transformation 
geomTransf Linear 1 
# Define Elements (columns and beam) 
# element elasticBeamColumn <element id> <start node> <end node> <area sq in> <E ksi> <I in4> <$transfTag> 
element elasticBeamColumn 1 1 2 14400.000000 3605.000000 108000000.000000 1 
element elasticBeamColumn 2 2 3 30720.000000 3605.000000 2684354560.000000 1 
element elasticBeamColumn 3 3 4 28672.000000 3605.000000 2505397589.000000 1 
element elasticBeamColumn 4 4 5 28672.000000 3605.000000 2505397589.000000 1 
element elasticBeamColumn 5 5 6 28672.000000 3605.000000 2505397589.000000 1 
element elasticBeamColumn 6 6 7 28672.000000 3605.000000 2505397589.000000 1 
