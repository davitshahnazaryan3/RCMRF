# Static analysis parameters
constraints Penalty 1.0e15 1.0e15;
numberer RCM;
system UmfPack;
test NormDispIncr 1e-08 6;
algorithm Newton;
integrator LoadControl 1.0;
analysis Static;
analyze 1;
loadConst -time 0.0;