wipe 
source mini_ida/model/model.tcl 
source mini_ida/model/loads.tcl 
source mini_ida/model/recorders.tcl 
wipeAnalysis 
constraints Transformation 
numberer RCM 
system Mumps 
test NormDispIncr 0.000001 100 
algorithm KrylovNewton 
integrator Newmark 0.5 0.25 
analysis Transient 
