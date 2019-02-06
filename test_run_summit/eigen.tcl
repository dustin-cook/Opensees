wipeAnalysis 
recorder Node -xml test_run_summit/mode_shape_1.xml -dT 1.000000 -node 6001  6004 -dof 1 3 "eigen 1" 
recorder Node -xml test_run_summit/mode_shape_2.xml -dT 1.000000 -node 6001  6004 -dof 1 3 "eigen 2" 
set numModes 4 
set lambda [eigen -fullGenLapack $numModes] 
set T {}
set pi [expr 2.0*asin(1.0)] 
foreach lam $lambda {
	lappend T [expr (2.0*$pi)/sqrt($lam)]
}
set period_file test_run_summit/period.txt 
set Periods [open $period_file "w"] 
foreach t $T { 
	puts $Periods " $t" 
} 
close $Periods 
integrator LoadControl 0 1 0 0 
test EnergyIncr 1.0e-10 100 0 
algorithm KrylovNewton 
numberer RCM 
constraints Transformation 
system ProfileSPD 
analysis Static 
set res [analyze 1] 
if {$res < 0} { 
    puts "Modal analysis failed" 
} 
remove recorders 
