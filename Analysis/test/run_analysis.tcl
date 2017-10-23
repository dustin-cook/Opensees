## Clear set up for this analysis 
wipe 
## Build Model and Analysis Parameters 
source Analysis/test/model.tcl 
source Analysis/test/recorders.tcl 
source Analysis/test/loads.tcl 
## ANALYSIS DEFINITION 
# Define Constraints 
constraints Plain 
# Define the DOF_numbered object 
numberer Plain 
# Construct Linear Solver and linear SOE Objects 
system BandGeneral 
# Construct Convergence Test 
test NormDispIncr 1.0e-6 6 
# Define Solution Algorithm 
algorithm Linear 
# Define Each Load Step (displacement controlled) 
integrator Newmark 0.5 0.25 
# Define analysis type 
analysis Transient 
## Run the Analysis 
analyze 3995 1.000000e-01 
puts "Done!" 
wipe 
