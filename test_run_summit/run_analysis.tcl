source outputs/simple_frame_3D/test/opensees_data/setup_analysis.tcl 
set singularity_check 0 
set collapse_check 0 
set currentStep [getTime] 
set ok 0 
set dt_reduce 1.000000 
set dt [expr 0.010000/$dt_reduce] 
puts "Analysis dt = $dt" 
set ok [analyze 2500 $dt] 
puts "analysis failure = $ok " 
if {$ok != 0} { 
puts "Analysis Failure: Convergence" 
wipe 
exit 
} 
if {$singularity_check == 1} { 
puts "Analysis Failure: Singularity" 
wipe 
exit 
} 
if {$collapse_check == 1} { 
puts "Analysis Failure: Collapse" 
wipe 
exit 
} 
if {$ok == 0} { 
puts "Analysis Complete!" 
set time_step_file outputs/simple_frame_3D/test/opensees_data/final_time_step_reduction.txt 
set TS [open $time_step_file "w"] 
puts $TS " $dt_reduce" 
close $TS 
wipe 
} 
