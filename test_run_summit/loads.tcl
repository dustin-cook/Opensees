pattern Plain 1 Linear {  
   load 6000 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6001 0.0 -150000.000000 0.0 0.0 0.0 0.0 
   load 6002 0.0 -150000.000000 0.0 0.0 0.0 0.0 
   load 6003 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6004 0.0 -150000.000000 0.0 0.0 0.0 0.0 
   load 6005 0.0 -150000.000000 0.0 0.0 0.0 0.0 
   load 6006 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6007 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6008 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6009 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6010 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6011 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6012 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6013 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6014 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6015 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6016 0.0 -0.000000 0.0 0.0 0.0 0.0 
   load 6017 0.0 -0.000000 0.0 0.0 0.0 0.0 
} 
constraints Transformation 
numberer RCM 
system BandGeneral 
test NormDispIncr 1.0e-5 1000 
algorithm Linear 
integrator LoadControl 0.1 
analysis Static 
set ok [analyze 10] 
analyze 10 
if {$ok != 0} { 
puts "Analysis Failure: Gravity Load Failure" 
wipe 
exit 
} 
puts "Gravity Load Complete" 
loadConst -time 0.0 
set dt 0.010000 
puts "EQ X dt = $dt" 
timeSeries Path 1 -dt $dt -filePath ground_motions/ICSB_1979/chan_13_accel_short.tcl -factor 386.000000 
pattern UniformExcitation 3 1 -accel 1 -fact 1.000000 
set lambda [eigen -fullGenLapack 3] 
set pi [expr 2.0*asin(1.0)] 
set i 0 
foreach lam $lambda {
    set i [expr $i+1] 
	set omega($i) [expr sqrt($lam)]
	set period($i) [expr 2*$pi/sqrt($lam)]
}
puts $period(1) 
puts $period(2) 
puts $period(3) 
set zeta 0.050000
set a0 [expr $zeta*2.0*$omega(1)*$omega(2)/($omega(1) + $omega(2))]
set a1 [expr $zeta*2.0/($omega(1) + $omega(2))]
rayleigh $a0 $a1 0.0 0.0 
