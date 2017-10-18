## Clear set up for this analysis 
wipe 
## Build Model and Analysis Parameters 
source Analysis/test/model.tcl 
source Analysis/test/recorders.tcl 
source Analysis/test/loads.tcl 
## ANALYSIS DEFINITION 
# Define Constraints 
constraints Transformation 
# Define the DOF_numbered object 
numberer RCM 
# Construct Linear Solver and linear SOE Objects 
system BandGeneral 
# Construct Convergence Test 
test NormDispIncr 1.0e-6 6 
# Define Solution ALgorithm 
algorithm Newton 
# Define Each Load Step (displacement controlled) 
integrator LoadControl 1 
# Define analysis type 
analysis Static 
## Run the Analysis 
analyze 1 
