# Define recorders 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_disp_x.txt -node 1  2  3  4  5  6  7 -dof 1 disp 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_disp_y.txt -node 1  2  3  4  5  6  7 -dof 2 disp 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_reaction_x.txt -node 1  2  3  4  5  6  7 -dof 1 reaction 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_reaction_y.txt -node 1  2  3  4  5  6  7 -dof 2 reaction 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_accel_x.txt -node 1  2  3  4  5  6  7 -dof 1 accel 
recorder Node -file outputs/ICBS_ns_wall_elastic/nodal_accel_y.txt -node 1  2  3  4  5  6  7 -dof 2 accel 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_1_force.txt -ele 1 force 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_2_force.txt -ele 2 force 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_3_force.txt -ele 3 force 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_4_force.txt -ele 4 force 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_5_force.txt -ele 5 force 
recorder Element -file outputs/ICBS_ns_wall_elastic/element_6_force.txt -ele 6 force 
# display displacement shape of the column 
recorder display "Displaced shape" 10 10 500 500 -wipe 
prp 200. 50. 1; 
vup 0 1 0; 
vpn 0 0 1; 
display 1 5 40 
