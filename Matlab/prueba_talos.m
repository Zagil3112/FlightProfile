% ...This program calculates the trajectory for a three stage sounding
% rocket.

%Declare global variables for use in every function

global g0 rho0 hscale Re Ma15C subs trans a b c d e f diam_talos A_talos Th_talos mdot_talos


deg			= pi/180;			% ...Convert degrees to radians

% ...Geophysical constants
g0			= 9.81;				% ...Sea level acceleration of gravity [m/s]
rho0			= 1.225;			% ...Sea level density of atmosphere
hscale			= 7.5e3;			% ...Density scale height [m]
Re			= 6378e3;			% ...Radius of the earth [m]

% Flow regime constants for the calculation of Drag
Ma15C			= 340.3;			% ...Speed of sound @ 15ºC [m/s]
subs			= 0.89;				% ...Limit of subsonic regime for this rocket nose (Hoerner 1965)
trans			= 1.13;				% ...High limit of transonic regime for this rocket nose (Hoerner 1965)

% Empirical coeficients for Drag calculation (from Hoerner, 1965)
a			= 0.1;		
b			= 1;
c			= 1.231;
d			= 0.141;
e			= 0.129;
f			= 1;

% ...Rocket constants for Talos-1 stage
diam_talos	= 0.76;					% ...Vehicle Diameter [m]
A_talos		= pi/4*(diam_talos)^2;			% ...Frontal area [m^2]
m0_talos	= 1996;					% ...Talos stage Fueled mass [kg]
mfinal_talos	= 496;					% ...Unfueled mass of Talos
Th_talos 	= 516e3;				% ...Thrust of Talos Stage [N]
mdot_talos	=-288.46;				% ...Rate of propellant consumption Talos stage [kg/s]
mprop_talos	= m0_talos - mfinal_talos;		% ...Talos Propellant mass [kg]
tburn_talos	= mprop_talos/abs(mdot_talos);		% ...Burn time Talos [s]


p = [g0 rho0 hscale Re Ma15C subs trans a b c d e f A_talos Th_talos mdot_talos];



t0 = 0;
tf = tburn_talos;
tspan = [t0 tf];

% ...Initial conditions:
v0		= 0;						% ...Initial velocity [m/s]
gamma0		= 85*deg;					% ...Initial flight path angle [rad]
x0		= 0;						% ...Initial downrange distance [km]
h0		= 0;						% ...Initial altitude
m0		= m0_talos;					% ...Lift off mass [kg]

z0=[v0 gamma0 x0 h0 m0];

options = odeset('Events',@ballis_event,'reltol',1e-12,'abstol',1e-12);	% ... SPECIFY TOLERANCE FOR THE ODE45 INTEGRATOR

[t,z]=ode45(@powered_talos,tspan,z0,options);

% Solution x(t) returned on the time interval [t0 tf] first flight phase:
v_1st		= z(:,1);	 				% ...Velocity Vector [Km/s]
gamma_1st	= z(:,2);	 				% ...Flight Path angle vector [degrees]
x_1st		= z(:,3);	 				% ...Downrange distance vector [Km]
h_1st		= z(:,4);	 				% ...Altitude vector [Km]
m_1st		= z(:,5);	 				% ...Mass variation [Kg]



% ...Initial conditions for the second flight phase (1st ballistic)

v0_2nd		= v_1st(end);					% ...Initial velocity (final 1st phase)
gamma0_2nd	= gamma_1st(end);				% ...Initial angle (final 1st phase)
x0_2nd		= x_1st(end);					% ...Initial downrange point (final 1st phase)
h0_2nd		= h_1st(end);					% ...Initial altitude (final 1st phase)
m0_2nd		= m_1st(end);					% ...Initial mass (Talos is detached)

t02	= tf;							% ...Initial time for 2nd phase = final time 1st phase
tf2	= t02+2000;		

tspan2 	= [t02, tf2];						%

z2_0 		= [v0_2nd gamma0_2nd x0_2nd h0_2nd m0_2nd];


[t2,z2,te,xe,ie]=ode45(@ballistic,tspan2,z2_0,options);

% Solution z2(t) returned on the time interval [t02 tf2]
v_2nd	   = z2(:,1);						% ...2nd phase (ballistic) velocity vector
gamma_2nd  = z2(:,2);						% ...2nd phase (ballistic) flight path angle vector
x_2nd	   = z2(:,3);						% ...2nd phase (ballistic) downrange vector
h_2nd	   = z2(:,4);						% ...2nd phase (ballistic) altitude vector
m_2nd	   = z2(:,5);						% ...Tomahawk + Nikha mass (costant)



plot(x_2nd,h_2nd)
