function Cd = Cd_reentry(mach)

mach_number = [0; 0.3; 0.5; 0.8; 0.81; 0.9; 0.96; 1.01; 1.2; 1.21; 1.6; 2; 3.5; 4.2; 6; 8; 99];
Cd_number = [.41 .41 .41 .45 0.49 .5 .58 .62 .71 0.71 .68 .66 .6587 .6482 .6275 .6 .6]';
Cd = interp1(mach_number, Cd_number, mach);
