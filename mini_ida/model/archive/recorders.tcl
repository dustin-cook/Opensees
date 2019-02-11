puts "Defining Recorders ..." 
recorder Node -xml mini_ida/model/nodal_disp_6000.xml -time -node 6000 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6000.xml -time -node 6000 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6001.xml -time -node 6001 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6001.xml -time -node 6001 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6002.xml -time -node 6002 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6002.xml -time -node 6002 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6003.xml -time -node 6003 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6003.xml -time -node 6003 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6004.xml -time -node 6004 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6004.xml -time -node 6004 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6005.xml -time -node 6005 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6005.xml -time -node 6005 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6006.xml -time -node 6006 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6006.xml -time -node 6006 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6007.xml -time -node 6007 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6007.xml -time -node 6007 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6008.xml -time -node 6008 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6008.xml -time -node 6008 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6009.xml -time -node 6009 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6009.xml -time -node 6009 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6010.xml -time -node 6010 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6010.xml -time -node 6010 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6011.xml -time -node 6011 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6011.xml -time -node 6011 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6012.xml -time -node 6012 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6012.xml -time -node 6012 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6013.xml -time -node 6013 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6013.xml -time -node 6013 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6014.xml -time -node 6014 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6014.xml -time -node 6014 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6015.xml -time -node 6015 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6015.xml -time -node 6015 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6016.xml -time -node 6016 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6016.xml -time -node 6016 -dof 1 3 accel 
recorder Node -xml mini_ida/model/nodal_disp_6017.xml -time -node 6017 -dof 1 3 disp 
recorder Node -xml mini_ida/model/nodal_accel_6017.xml -time -node 6017 -dof 1 3 accel
recorder display "Displaced shape" 10 10 500 500 -wipe 
prp 200.0 50.0 50.0; 
vup 0.0 1.0 0.0; 
vpn 0.4 0.25 1; 
display 1 5 1.000000 
puts "Defining Recorders Complete" 
