# Define Gravity Loads (node id, axial, shear, moment) 
timeSeries Linear 1 
pattern Plain 1 1 {  
   load 1 0.0 -0.000000 0.0 
   load 2 0.0 -2200.000000 0.0 
   load 3 0.0 -2325.000000 0.0 
   load 4 0.0 -2325.000000 0.0 
   load 5 0.0 -2325.000000 0.0 
   load 6 0.0 -2325.000000 0.0 
   load 7 0.0 -1900.000000 0.0 
} 
## GRAVITY ANALYSIS 
constraints Plain 
numberer Plain 
system BandGeneral 
integrator LoadControl 0.1 
analysis Static	 
analyze 10 
loadConst -time 0.0 
# Define Load Pattern 
pattern Plain 2 Linear { 
  load 1 0.000000 0.0 0.0 
  load 2 0.000000 0.0 0.0 
  load 3 0.000000 0.0 0.0 
  load 4 0.000000 0.0 0.0 
  load 5 0.000000 0.0 0.0 
  load 6 0.000000 0.0 0.0 
  load 7 0.000000 0.0 0.0 
} 
# Define Seismic Excitation Pattern 
timeSeries Path 2 -dt 0.010000 -filePath ground_motions/ICSB_1979/eq_ew_ground.tcl -factor 386. 
pattern UniformExcitation 3 1 -accel 2 
# set damping based on first eigen mode 
set lambda [expr [eigen -fullGenLapack 1]] 
set omega [expr sqrt($lambda)] 
rayleigh 0 0 0 [expr 2*5.000000e-02/$omega] 
rayleigh 0 0 0 [expr 2*5.000000e-02/$omega] 
