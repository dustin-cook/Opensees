puts "Defining Ground Motions and Scale ..." 
timeSeries Path 1 -dt $dt -filePath ground_motions/ICSB_1979/chan_13_accel_15.tcl -factor 386.000000 
pattern UniformExcitation 3 1 -accel 1 -fact 3.000000 
timeSeries Path 2 -dt $dt -filePath ground_motions/ICSB_1979/chan_11_accel_15.tcl -factor 386.000000 
pattern UniformExcitation 4 3 -accel 2 -fact 3.000000 
