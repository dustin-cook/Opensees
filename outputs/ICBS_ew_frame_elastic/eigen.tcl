# Eigen analysis - 10.4 from "Dynamics of Structures" book by Anil Chopra - using equalDOF and very high Ib 
wipeAnalysis 
# record eigenvectors 
recorder Node -file Aoutputs/ICBS_ew_frame_elastic/mode_shape_1.txt -dT 0.020000 -node 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42 -dof 1 "eigen 1" 
# perform eigen analysis
set numModes 6 
set lambda [eigen -fullGenLapack $numModes] 
set T {}
set pi 3.141593
foreach lam $lambda {
	lappend T [expr (2*$pi)/sqrt($lam)]
}
# write the output file cosisting of periods 
set period_file outputs/ICBS_ew_frame_elastic/period.txt 
set Periods [open $period_file "w"] 
foreach t $T { 
	puts $Periods " $t" 
} 
close $Periods 
# Run a one step gravity load with no loading (to record eigenvectors)
integrator LoadControl 0 1 0 0 
test EnergyIncr 1.0e-10 100 0 
algorithm Newton 
numberer RCM 
constraints Transformation 
system ProfileSPD 
analysis Static 
set res [analyze 1] 
if {$res < 0} { 
    puts "Modal analysis failed" 
} 
remove recorders 
