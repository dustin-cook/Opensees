# Define recorders 
recorder Node -file Analysis/test/nodal_disp_x.txt -node 1 2 3 4 -dof 1 disp 
recorder Node -file Analysis/test/nodal_disp_y.txt -node 1 2 3 4 -dof 2 disp 
recorder Node -file Analysis/test/nodal_reaction_x.txt -node 1 2 3 4 -dof 1 reaction 
recorder Node -file Analysis/test/nodal_reaction_y.txt -node 1 2 3 4 -dof 2 reaction 
recorder Element -file Analysis/test/element_1_force.txt -ele 1 force 
recorder Element -file Analysis/test/element_2_force.txt -ele 2 force 
recorder Element -file Analysis/test/element_3_force.txt -ele 3 force 
