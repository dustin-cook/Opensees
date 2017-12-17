## Clear set up for this analysis 
wipe 
## Build Model and Analysis Parameters 
source outputs/ICBS_ew_frame_elastic/model.tcl 
source outputs/ICBS_ew_frame_elastic/eigen.tcl 
source outputs/ICBS_ew_frame_elastic/loads.tcl 
source outputs/ICBS_ew_frame_elastic/recorders.tcl 
## ANALYSIS DEFINITION 
wipeAnalysis 
# Define Constraints 
constraints Plain 
# Define the DOF_numbered object 
numberer Plain 
# Construct Linear Solver and linear SOE Objects 
system BandGeneral 
# Define Solution Algorithm 
algorithm Linear 
# Define Each Load Step (displacement controlled) 
integrator Newmark 0.5 0.25 
# Define analysis type 
analysis Transient 
## Run the Analysis 
analyze 5688 0.010000 
puts "Done!" 
wipe 
