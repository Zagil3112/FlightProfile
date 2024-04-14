function rho = atmosphere1(h);

rho0 = 1.225;
hscale = 7.5e3;

rho = rho0*exp(-h/hscale);
