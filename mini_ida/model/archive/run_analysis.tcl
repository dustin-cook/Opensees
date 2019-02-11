source mini_ida/model/setup_analysis.tcl 
set singularity_check 0 
set collapse_check 0 
set currentStep [getTime] 
set ok 0 
set converge_tol_file [open mini_ida/model/converge_tol_file.txt w] 
set converge_file [open mini_ida/model/converge_file.txt w] 
while {$ok == 0 && $currentStep < 15.000000 && $collapse_check == 0 && $singularity_check == 0} { 
puts "Progress: $currentStep out of 15.000000" 
set tol 0.000010 
test NormDispIncr $tol 10 
algorithm KrylovNewton 
set dt_reduce 1.000000 
set dt [expr 0.010000/$dt_reduce] 
set dt_max 0.010000 
set dt_min [expr 0.010000/($dt_reduce*100)] 
set ok [analyze 1 $dt $dt_min $dt_max] 
puts "analysis failure = $ok " 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000010, step_length/1.000000" 
set tol 0.000010 
test NormDispIncr $tol 10 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000010, step_length/10.000000" 
set tol 0.000010 
test NormDispIncr $tol 10 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000010, step_length/20.000000" 
set tol 0.000010 
test NormDispIncr $tol 10 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000010, step_length/50.000000" 
set tol 0.000010 
test NormDispIncr $tol 10 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000100, step_length/1.000000" 
set tol 0.000100 
test NormDispIncr $tol 10 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000100, step_length/10.000000" 
set tol 0.000100 
test NormDispIncr $tol 10 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000100, step_length/20.000000" 
set tol 0.000100 
test NormDispIncr $tol 10 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.000100, step_length/50.000000" 
set tol 0.000100 
test NormDispIncr $tol 10 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.001000, step_length/1.000000" 
set tol 0.001000 
test NormDispIncr $tol 10 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.001000, step_length/10.000000" 
set tol 0.001000 
test NormDispIncr $tol 10 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.001000, step_length/20.000000" 
set tol 0.001000 
test NormDispIncr $tol 10 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.001000, step_length/50.000000" 
set tol 0.001000 
test NormDispIncr $tol 10 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.010000, step_length/1.000000" 
set tol 0.010000 
test NormDispIncr $tol 10 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.010000, step_length/10.000000" 
set tol 0.010000 
test NormDispIncr $tol 10 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.010000, step_length/20.000000" 
set tol 0.010000 
test NormDispIncr $tol 10 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.010000, step_length/50.000000" 
set tol 0.010000 
test NormDispIncr $tol 10 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.100000, step_length/1.000000" 
set tol 0.100000 
test NormDispIncr $tol 1000 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.100000, step_length/10.000000" 
set tol 0.100000 
test NormDispIncr $tol 1000 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.100000, step_length/20.000000" 
set tol 0.100000 
test NormDispIncr $tol 1000 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 0.100000, step_length/50.000000" 
set tol 0.100000 
test NormDispIncr $tol 1000 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 1.000000, step_length/1.000000" 
set tol 1.000000 
test NormDispIncr $tol 1000 
set step_reduce 1.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 1.000000, step_length/10.000000" 
set tol 1.000000 
test NormDispIncr $tol 1000 
set step_reduce 10.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 1.000000, step_length/20.000000" 
set tol 1.000000 
test NormDispIncr $tol 1000 
set step_reduce 20.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
if {$ok != 0} { 
puts "analysis failed, try try tolerance = 1.000000, step_length/50.000000" 
set tol 1.000000 
test NormDispIncr $tol 1000 
set step_reduce 50.000000 
set dt [expr 0.010000/$step_reduce] 
set ok [analyze 1 $dt] 
} 
set currentStep [getTime] 
set converge_tol_log "$currentStep $tol" 
puts $converge_tol_file $converge_tol_log 
if {$ok == 0} { 
set node_at_floor_1 6001 
set floor_displ_1 "[nodeDisp $node_at_floor_1 1]" 
puts "First Story Disp = $floor_displ_1" 
set height_floor_1 120.000000 
set floor_drift_1 [expr abs($floor_displ_1/$height_floor_1)] 
set check_QNAN_1 [string first QNAN $floor_displ_1 1] 
set check_IND_1 [string first IND $floor_displ_1 1] 
if {($floor_displ_1 > 1000000) || ($check_QNAN_1 != -1) || ($check_IND_1 != -1)} { 
set singularity_check 1 
} 
if {$floor_drift_1 > 0.060000} { 
set collapse_check 1 
} 
} 
} 
if {$currentStep > 14.700000 || $singularity_check == 1 || $collapse_check == 1} { 
puts "Analysis Fully Converged" 
puts $converge_file 1 
} else { 
puts "Analysis NOT FULLY Converged" 
puts $converge_file 0 
} 
close $converge_tol_file 
close $converge_file 
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
set time_step_file mini_ida/model/final_time_step_reduction.txt 
set TS [open $time_step_file "w"] 
puts $TS " $dt_reduce" 
close $TS 
wipe 
} 
