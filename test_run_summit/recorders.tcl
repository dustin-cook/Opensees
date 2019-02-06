recorder Node -xml outputs/simple_frame_3D/test/opensees_data/nodal_disp_x.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 1 disp 
recorder Node -xml outputs/simple_frame_3D/test/opensees_data/nodal_accel_x.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 1 accel 
recorder Node -xml outputs/simple_frame_3D/test/opensees_data/nodal_disp_z.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 3 disp 
recorder Node -xml outputs/simple_frame_3D/test/opensees_data/nodal_accel_z.xml -time -node 6000  6001  6002  6003  6004  6005  6006  6007  6008  6009  6010  6011  6012  6013  6014  6015  6016  6017 -dof 3 accel 
recorder Element -xml outputs/simple_frame_3D/test/opensees_data/element_force.xml -time -ele 1  2  3  4  5  6 -dof 1 2 3 6 12 localForce 
recorder display "Displaced shape" 10 10 500 500 -wipe 
prp 200.0 50.0 50.0; 
vup 0.0 1.0 0.0; 
vpn 0.4 0.25 1; 
display 1 5 1.000000 
