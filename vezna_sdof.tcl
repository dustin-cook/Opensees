wipe; 
model basic -ndm 2 -ndf 3; 
file mkdir data; 
node 1 0. 0.; 
node 2 0. 432. 
fix 1 1 1 1; 
mass 2 5.18 0. 0.; 
geomTransf Linear 1; 
element elasticBeamColumn 1 1 2 3600 3225 4792.896333 1; 
recorder Node -file vezna_accel.txt -time -node 2 -dof 1 accel; 
timeSeries Linear 1 
pattern Plain 1 1 { 
   load 2 0. -2000. 0.; 
} 
constraints Plain; 
numberer Plain; 
system BandGeneral; 
algorithm Linear; 
integrator LoadControl 0.1; 
analysis Static 
analyze 10; 
loadConst -time 0.0; 
set G 386 
timeSeries Path 2 -dt 0.005000 -filePath A10000.tcl -factor $G; 
pattern UniformExcitation 2 1 -accel 2; 
set freq [expr [eigen -fullGenLapack 1]**0.5] 
set dampRatio 0.02 
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq] 
wipeAnalysis; 
constraints Plain; 
numberer Plain; 
system BandGeneral; 
algorithm Linear 
integrator Newmark 0.5 0.25 ; 
analysis Transient; 
analyze 3995 0.010000; 
puts "Done!" 
wipe 
