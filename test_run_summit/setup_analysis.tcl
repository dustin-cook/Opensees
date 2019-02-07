wipe 
source test_run_summit/model.tcl 
source test_run_summit/eigen.tcl 
source test_run_summit/loads.tcl 
source test_run_summit/recorders.tcl 
wipeAnalysis 
constraints Transformation 
numberer RCM 
system BandGeneral 
test NormDispIncr 0.000001 100 
algorithm KrylovNewton 
integrator Newmark 0.5 0.25 
analysis Transient 
