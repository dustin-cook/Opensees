wipe 
source outputs/simple_frame_3D/test/opensees_data/model.tcl 
source outputs/simple_frame_3D/test/opensees_data/eigen.tcl 
source outputs/simple_frame_3D/test/opensees_data/loads.tcl 
source outputs/simple_frame_3D/test/opensees_data/recorders.tcl 
wipeAnalysis 
constraints Transformation 
numberer RCM 
system BandGeneral 
test NormDispIncr 0.000001 100 
algorithm KrylovNewton 
integrator Newmark 0.5 0.25 
analysis Transient 
