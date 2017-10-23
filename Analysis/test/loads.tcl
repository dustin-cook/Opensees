# Define Gravity Loads (node id, axial, shear, moment) 
pattern Plain 1 Linear {  
   load 3 0.0 0 0.0 
   load 4 0.0 0 0.0 
} 
# Define Load Pattern 
pattern Plain 2 Linear { 
  load 3 0 0.0 0.0 
} 
# Define Seismic Excitation Pattern 
timeSeries Path 1 -dt 0.005 -filePath ground_motions/A10000.tcl -factor 386 
pattern UniformExcitation 3 1 -accel 1 
# set damping based on first eigen mode 
set freq [expr [eigen -fullGenLapack 1]**0.5] 
  rayleigh 0. 0. 0. [expr 2*5.000000e-02/$freq] 
