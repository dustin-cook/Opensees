recorder Node -xml test_run_summit/nodal_disp_x.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 1 disp 
recorder Node -xml test_run_summit/nodal_accel_x.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 1 accel 
recorder Node -xml test_run_summit/nodal_disp_z.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 3 disp 
recorder Node -xml test_run_summit/nodal_accel_z.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 3 accel 
recorder Element -xml test_run_summit/element_force.xml -time -ele 1  2  3  4  5  6 -dof 1 2 3 6 12 localForce 
