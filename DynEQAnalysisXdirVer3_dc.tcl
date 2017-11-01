# # Dynamic Analysis

puts "EQ Analysis start"

# source ReadSMDfile.tcl
source Setup.tcl;
# source Accel.tcl

set GMdirection1 1;
# set dt 0.005; # change
# set GMfact 8.1467;

# set DtAnalysis [expr $dt/5.];
# set TmaxAnalysis 40.

# # Damping
set xDamp 0.02;
set MpropSwitch 0.0;
set KcurrSwitch 0.0;
set KcommSwitch 0.8;
set KinitSwitch 0.55;
set omegaI [expr [eigen -fullGenLapack 1]**0.5]
set alphaM [expr $MpropSwitch*$xDamp*(2.)];	# M-prop. damping; D = alphaM*M
set betaKcurr [expr $KcurrSwitch*2.*$xDamp/($omegaI)];         		# current-K;      +beatKcurr*KCurrent
set betaKcomm [expr $KcommSwitch*2.*$xDamp/($omegaI)];   		# last-committed K;   +betaKcomm*KlastCommitt
set betaKinit [expr $KinitSwitch*2.*$xDamp/($omegaI)];         			# initial-K;     +beatKinit*Kini

# region 1 -eleRange 12 119 -rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm ; 				# RAYLEIGH damping
rayleigh $alphaM $betaKcurr $betaKinit $betaKcomm ; 				# RAYLEIGH damping

puts "Damping applied"

# # Analysis Parameters

set constraintsTypeDynamic Transformation;
constraints $constraintsTypeDynamic ;

set numbererTypeDynamic RCM
numberer $numbererTypeDynamic 

set systemTypeDynamic BandGeneral;	# try UmfPack for large problems
system $systemTypeDynamic

set TolDynamic 1.0e-8;
set Tol 1.0e-8;
set Tol_adv 1.0e-5;                        # Convergence Test: tolerance
set maxNumIterDynamic 100;                # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
set printFlagDynamic 0;                # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
set testTypeDynamic EnergyIncr;	# Convergence-test type EnergyIncr
test $testTypeDynamic $TolDynamic $maxNumIterDynamic $printFlagDynamic;
# for improved-convergence procedure:
	set maxNumIterConvergeDynamic 1500;	
	set printFlagConvergeDynamic 0;
	
set algorithmTypeDynamic Newton 
algorithm $algorithmTypeDynamic;

set NewmarkGamma 0.5;	# Newmark-integrator gamma parameter 
set NewmarkBeta 0.25;	# Newmark-integrator beta parameter
set integratorTypeDynamic Newmark;
integrator $integratorTypeDynamic $NewmarkGamma $NewmarkBeta -maxDU 0.01
# integrator TRBDF2

set analysisTypeDynamic Transient
analysis $analysisTypeDynamic 

#  ---------------------------------    perform Dynamic Ground-Motion Analysis
# the following commands are unique to the Uniform Earthquake excitation
set IDloadTag1 400;	# for uniformSupport excitation
# Uniform EXCITATION: acceleration input

puts "Starts dynamic analysis"
set GMfatt 386;		# data in input file is in cm Unifts -- ACCELERATION TH
set AccelSeries1 "Series -dt $dt -filePath A10000.tcl -factor  $GMfatt";	# time series information
pattern UniformExcitation  $IDloadTag1  $GMdirection1 -accel  $AccelSeries1  ;		# create Unifform excitation
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];
set ok [analyze $Nsteps $DtAnalysis];	
puts "Analysis $ok"

if {$ok != 0} {      ;					# analysis was not successful.
	# --------------------------------------------------------------------------------------------------
	# change some analysis parameters to achieve convergence
	# performance is slower inside this loop
	#    Time-controlled analysis
	puts "Advanced analysis"
	set ok 0;
	set controlTime [getTime];
	while {$controlTime < $TmaxAnalysis && $ok == 0} {
		set controlTime [getTime]
		set ok [analyze 1 $DtAnalysis]
		if {$ok != 0} {
			puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Change tolerance
		if {$ok != 0} {
		puts "Change tolerance at normal time step"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 1 $DtAnalysis]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 2
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/2.]
		puts "Time Step divided by 2"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 2 and tolarence
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/2.]
		puts "Time Step divided by 2"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 2 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 5
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/5.]
		puts "Time Step divided by 5"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 5 and tolerance
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/5.]
		puts "Time Step divided by 5 and tol"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 5 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 10
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/10.]
		puts "Time Step divided by 10"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 10 and tol
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/10.]
		puts "Time Step divided by 10 and tol"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 10 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 100
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/100.]
		puts "Time Step divided by 100"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 100 and tol
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/100.]
		puts "Time Step divided by 100"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 1000
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/1000.]
		puts "Time Step divided by 1000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 1000 and tolerance
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/1000.]
		puts "Time Step divided by 1000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
        # Changigng time step by 10000
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/10000.]
		puts "Time Step divided by 10000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 10000 and tol
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/10000.]
		puts "Time Step divided by 10000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		
		# Changigng time step by 100000
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/100000.]
		puts "Time Step divided by 100000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}	
		
		# Changigng time step by 100000 and tol
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/100000.]
		puts "Time Step divided by 100000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}	

		# Changigng time step by 1000000
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/1000000.]
		puts "Time Step divided by 1000000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}		
		
		# Changigng time step by 1000000 and tol
		if {$ok != 0} {
		set DtCurrent [expr $DtAnalysis/1000000.]
		puts "Time Step divided by 1000000"
		puts "Trying Modified Newton .."
			test EnergyIncr   $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic	
		}
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm NewtonLineSearch -tol 0.8 -maxiter $maxNumIterConvergeDynamic
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
		puts "Trying Newton with Initial Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initial
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}
		if {$ok != 0} {
			puts "Trying Newton with Initial-then-Current Tangent .."
			test EnergyIncr  $Tol_adv $maxNumIterConvergeDynamic  0
			algorithm Newton -initialThenCurrent
			set ok [analyze 100 $DtCurrent]
			test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			algorithm $algorithmTypeDynamic
		}		
		
		
		
		
		# Changigng time step by 1000
		# if {$ok != 0} {
		# set DtCurrent [expr $DtAnalysis/1000.]
		# puts "Time Step divided by 1000"
		# puts "Trying Modified Newton .."
			# test EnergyIncr   $Tol_adv 200  0
			# algorithm Newton
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic	
		# }
		# if {$ok != 0} {
		# puts "Trying Newton with Initial Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initial
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		# if {$ok != 0} {
			# puts "Trying Newton with Initial-then-Current Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initialThenCurrent
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		
        # Changigng time step by 10000
		# if {$ok != 0} {
		# set DtCurrent [expr $DtAnalysis/10000.]
		# puts "Time Step divided by 10000"
		# puts "Trying Modified Newton .."
			# test EnergyIncr   $Tol_adv 200  0
			# algorithm Newton
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic	
		# }
		# if {$ok != 0} {
		# puts "Trying Newton with Initial Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initial
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		# if {$ok != 0} {
			# puts "Trying Newton with Initial-then-Current Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initialThenCurrent
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		
		# Changigng time step by 100000
		# if {$ok != 0} {
		# set DtCurrent [expr $DtAnalysis/100000.]
		# puts "Time Step divided by 100000"
		# puts "Trying Modified Newton .."
			# test EnergyIncr   $Tol_adv 200  0
			# algorithm Newton
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic	
		# }
		# if {$ok != 0} {
		# puts "Trying Newton with Initial Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initial
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		# if {$ok != 0} {
			# puts "Trying Newton with Initial-then-Current Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initialThenCurrent
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }	

		# Changigng time step by 1000000
		# if {$ok != 0} {
		# set DtCurrent [expr $DtAnalysis/1000000.]
		# puts "Time Step divided by 1000000"
		# puts "Trying Modified Newton .."
			# test EnergyIncr   $Tol_adv 200  0
			# algorithm Newton
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic	
		# }
		# if {$ok != 0} {
		# puts "Trying Newton with Initial Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initial
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }
		# if {$ok != 0} {
			# puts "Trying Newton with Initial-then-Current Tangent .."
			# test EnergyIncr  $Tol_adv 200  0
			# algorithm Newton -initialThenCurrent
			# set ok [analyze 100 $DtCurrent]
			# test $testTypeDynamic $TolDynamic $maxNumIterDynamic  0
			# algorithm $algorithmTypeDynamic
		# }		
			
	}
};      # end if ok !0
set An "Analysis.txt"
set Anal [open $An "w"]
puts $Anal "$ok" 
close $Anal
# puts "Ground Motion Done. End Time: [getTime]"