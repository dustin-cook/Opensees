# Define Gravity Loads (node id, axial, shear, moment) 
timeSeries Linear 1 
pattern Plain 1 1 {  
   load 1 0.0 -0.000000 0.0 
   load 2 0.0 -0.000000 0.0 
   load 3 0.0 -0.000000 0.0 
   load 4 0.0 -0.000000 0.0 
   load 5 0.0 -0.000000 0.0 
   load 6 0.0 -0.000000 0.0 
   load 7 0.0 -220.000000 0.0 
   load 8 0.0 -440.000000 0.0 
   load 9 0.0 -440.000000 0.0 
   load 10 0.0 -440.000000 0.0 
   load 11 0.0 -440.000000 0.0 
   load 12 0.0 -220.000000 0.0 
   load 13 0.0 -232.500000 0.0 
   load 14 0.0 -465.000000 0.0 
   load 15 0.0 -465.000000 0.0 
   load 16 0.0 -465.000000 0.0 
   load 17 0.0 -465.000000 0.0 
   load 18 0.0 -232.500000 0.0 
   load 19 0.0 -232.500000 0.0 
   load 20 0.0 -465.000000 0.0 
   load 21 0.0 -465.000000 0.0 
   load 22 0.0 -465.000000 0.0 
   load 23 0.0 -465.000000 0.0 
   load 24 0.0 -232.500000 0.0 
   load 25 0.0 -232.500000 0.0 
   load 26 0.0 -465.000000 0.0 
   load 27 0.0 -465.000000 0.0 
   load 28 0.0 -465.000000 0.0 
   load 29 0.0 -465.000000 0.0 
   load 30 0.0 -232.500000 0.0 
   load 31 0.0 -232.500000 0.0 
   load 32 0.0 -465.000000 0.0 
   load 33 0.0 -465.000000 0.0 
   load 34 0.0 -465.000000 0.0 
   load 35 0.0 -465.000000 0.0 
   load 36 0.0 -232.500000 0.0 
   load 37 0.0 -190.000000 0.0 
   load 38 0.0 -380.000000 0.0 
   load 39 0.0 -380.000000 0.0 
   load 40 0.0 -380.000000 0.0 
   load 41 0.0 -380.000000 0.0 
   load 42 0.0 -190.000000 0.0 
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
  load 8 0.000000 0.0 0.0 
  load 9 0.000000 0.0 0.0 
  load 10 0.000000 0.0 0.0 
  load 11 0.000000 0.0 0.0 
  load 12 0.000000 0.0 0.0 
  load 13 0.000000 0.0 0.0 
  load 14 0.000000 0.0 0.0 
  load 15 0.000000 0.0 0.0 
  load 16 0.000000 0.0 0.0 
  load 17 0.000000 0.0 0.0 
  load 18 0.000000 0.0 0.0 
  load 19 0.000000 0.0 0.0 
  load 20 0.000000 0.0 0.0 
  load 21 0.000000 0.0 0.0 
  load 22 0.000000 0.0 0.0 
  load 23 0.000000 0.0 0.0 
  load 24 0.000000 0.0 0.0 
  load 25 0.000000 0.0 0.0 
  load 26 0.000000 0.0 0.0 
  load 27 0.000000 0.0 0.0 
  load 28 0.000000 0.0 0.0 
  load 29 0.000000 0.0 0.0 
  load 30 0.000000 0.0 0.0 
  load 31 0.000000 0.0 0.0 
  load 32 0.000000 0.0 0.0 
  load 33 0.000000 0.0 0.0 
  load 34 0.000000 0.0 0.0 
  load 35 0.000000 0.0 0.0 
  load 36 0.000000 0.0 0.0 
  load 37 0.000000 0.0 0.0 
  load 38 0.000000 0.0 0.0 
  load 39 0.000000 0.0 0.0 
  load 40 0.000000 0.0 0.0 
  load 41 0.000000 0.0 0.0 
  load 42 0.000000 0.0 0.0 
} 
# Define Seismic Excitation Pattern 
timeSeries Path 2 -dt 0.010000 -filePath ground_motions/ICSB_1979/eq_ew_ground.tcl -factor 386. 
pattern UniformExcitation 3 1 -accel 2 
# set damping based on first eigen mode 
set lambda [expr [eigen -fullGenLapack 1]] 
set omega [expr sqrt($lambda)] 
rayleigh 0 0 0 [expr 2*5.000000e-02/$omega] 
rayleigh 0 0 0 [expr 2*5.000000e-02/$omega] 
