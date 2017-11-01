# Define Gravity Loads (node id, axial, shear, moment) 
timeSeries Linear 1 
pattern Plain 1 1 {  
   load 1 0.0 -0.000000 0.0 
   load 2 0.0 -386.000000 0.0 
} 
## GRAVITY ANALYSIS 
constraints Plain 
numberer Plain 
system BandGeneral 
integrator LoadControl 0.1 
analysis Static	 
analyze 10 
loadConst -time 0.0

