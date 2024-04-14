% ...THIS PROGRAM CALCULATES THE TRAJECTORY FOR A THREE STAGE SOUNDING ROCKET

% Declare global variables for use in every function

global rho0 Re R_polar subs trans a b c d e f omega diam_talos A_talos Th_talos mdot_talos diam_toma A_toma Th_toma mdot_toma diam_nikha A_nikha Th_nikha mdot_nikha epsilon mu J2 J3 J4

deg	= pi/180;			% Convert degrees to radians
rho0	= 1.225;			% Sea level density of atmosphere

% NGA/NASA EGM96, N=M=360, WGS84 REVISED EARTH GRAVITATIONAL MODEL
Re	= 6378.137e3;			% Semi-major axis of WGS 84 Ellipsoid (Radius of the earth equator for the NGA/NASA EGM96 WGS84) 	   						% Gravitational model[m]
R_polar	= 6356.752314245e3;		% Semi-minor axis of WGS 84 Ellipsoid (Polar radius of the earth) [m]
epsilon	= 1.0/298.2572235630;		% Flattening of WGS 84 Ellipsoid
omega	= 7.2921150e-5;			% Angular speed of Earth's rotation [rad/sec]
mu	= 3.986004418e14; 		% mu=GMe
J2	= 1.08263e-3;			% J2, J3, J4 harmonics for modelling the WGS84
J3	= 2.532153e-7;
J4	= 1.6109876e-7;

% FLOW REGIME CONSTANTS FOR DRAG CALCULATION

subs	= 0.89;				% Limit of subsonic regime for this rocket nose (Hoerner 1965)
trans	= 1.13;				% High limit of transonic regime for this rocket nose (Hoerner 1965)

% Empirical coeficients for Drag calculation (from Hoerner, 1965)
a	= 0.1;		
b	= 1;
c	= 1.231;
d	= 0.141;
e	= 0.129;
f	= 1;

% ROCKET CONSTANTS FOR TALOS-1 STAGE
diam_talos	= 0.76;				% Vehicle Diameter [m]
A_talos		= pi/4*(diam_talos)^2;		% Frontal area [m^2]
m0_talos	= 1996;				% Talos stage Fueled mass [kg]
mfinal_talos	= 496;				% Talos Unfueled mass
Th_talos 	= 516e3;			% Thrust of Talos Stage [N]
mdot_talos	=-288.46;			% Rate of propellant consumption Talos stage [kg/s]
mprop_talos	= m0_talos - mfinal_talos;	% Talos Propellant mass [kg]
tburn_talos	= mprop_talos/abs(mdot_talos);	% Burn time Talos [s]

% ROCKET CONSTANTS FOR TOMAHAWK STAGE
diam_toma	= 0.58;				% Vehicle Diameter 
A_toma		= pi/4*(diam_toma)^2;		% Frontal area [m^2]
m0_toma		= 1363;				% Tomahawk fueled mass [kg]
mfinal_toma	= 602;				% Tomahawk Unfueled mass
Th_toma		= 457e3;			% Thrust of Tomahawk stage [N]
mdot_toma	= -217.43;			% Toma Rate of propellant comsuption [kg/s]
mprop_toma	= m0_toma - mfinal_toma;	% Tomahaw propellant mass [kg]
tburn_toma	= mprop_toma/abs(mdot_toma);	% Burn Time Tomahawk [s]

% ROCKET CONSTANTS FOR NIKHA STAGE
diam_nikha	= 0.44;				% Vehicle Diameter [m]
A_nikha		= pi/4*(diam_nikha)^2;		% Frontal area [m^2]
m0_nikha	= 399;				% Nikha fueled mass [kg]
mfinal_nikha	= 70;				% Nikha unfueled mass [kg]
Th_nikha	= 50500;			% Thrust of Nikha stage [N]
mdot_nikha	= -19.35;			% Nikha propellant consumption [Kg/s]
mprop_nikha	= m0_nikha - mfinal_nikha;	% Nikha propellant mass [kg]
tburn_nikha	= mprop_nikha/abs(mdot_nikha);	% Burn time Nikha

% PAYLOAD MASS
m_pay		= 120; % [kg]

% INITIAL CONDITIONS
t0 = 0;
tf1 = tburn_talos;
tspan1 = [t0 tf1];

v1st_0		= 0;					% Initial velocity
A1st_0		= 90*deg;				% Initial azimuth
phi1st_0	= 85*deg;				% Initial angle wrt local horizontal
delta1st_0	= 37.940194*deg;			% Initial latitude
lamda1st_0	= -75.466389*deg;			% Initial longitude
r1st_0		= Re*(1-epsilon*(sin(delta1st_0))^2);	% Initial datum sea level due to WGS84 model
m1st_0		= m0_talos+m0_toma+m0_nikha+m_pay;	% Total initial mass

z1_0=[v1st_0 A1st_0 phi1st_0 r1st_0 delta1st_0 lamda1st_0 m1st_0];

options1 = odeset('Events',@ballis_event,'reltol',1e-12,'abstol',1e-12);  % Specify tolerance and event stop for ode45 integrator
options2 = odeset('Events',@reentry_event,'reltol',1e-12,'abstol',1e-12);

[t1,z1]=ode45(@powered_talos,tspan1,z1_0,options1);

% Solution x(t) returned on the time interval [t0 tf] first flight phase:
v1st		= z1(:,1);	 % Velocity [m/s]
A1st		= z1(:,2);	 % Azimuth angle [rad]
phi1st		= z1(:,3);	 % Flight path angle [rad]
r1st		= z1(:,4);	 % Distance from Earth center [m]
delta1st	= z1(:,5);	 % Latitude [rad]
lamda1st	= z1(:,6);	 % Longitude [rad]
m1st		= z1(:,7);	 % Mass [kg]

% Initial conditions for the second flight phase (1st ballistic)

v2nd_0		= v1st(end);				
A2nd_0		= A1st(end);				
phi2nd_0	= phi1st(end);					
r2nd_0		= r1st(end);					
delta2nd_0	= delta1st(end);					
lamda2nd_0	= lamda1st(end);
m2nd_0		= m0_toma+m0_nikha+m_pay; 	% Talos stage has been detached

t02	= tf1;		% Initial time for 2nd phase = final time 1st phase
tf2	= t02+2;	% Two seconds between the stages	
tspan2 	= [t02, tf2];						%

z2_0 		= [v2nd_0 A2nd_0 phi2nd_0 r2nd_0 delta2nd_0 lamda2nd_0 m2nd_0];

[t2,z2]=ode45(@ballistic,tspan2,z2_0,options1);

% Solution z2(t) returned on the time interval [t02 tf2]
v2nd	   	= z2(:,1);	% Velocity [m/s]
A2nd  		= z2(:,2);	% Azimuth angle [rad]
phi2nd	   	= z2(:,3);	% Flight path angle
r2nd	   	= z2(:,4);	% Distance from Earth Center
delta2nd	= z2(:,5);	% Latitude
lamda2nd	= z2(:,6);	% Longitude
m2nd		= z2(:,7);	% Mass

% Initial conditions for the third flight phase (powered)

v3rd_0		= v2nd(end);					
A3rd_0		= A2nd(end);				
phi3rd_0	= phi2nd(end);					
r3rd_0		= r2nd(end);					
delta3rd_0	= delta2nd(end);					
lamda3rd_0	= lamda2nd(end);
m3rd_0		= m2nd(end);

t03	= tf2;			% Initial time for 3rd phase = final time 2nd phase
tf3	= t03+tburn_toma;
tspan3	= [t03 tf3];

z3_0		= [v3rd_0 A3rd_0 phi3rd_0 r3rd_0 delta3rd_0 lamda3rd_0 m3rd_0];

[t3,z3]=ode45(@powered_tomahawk,tspan3,z3_0,options1);

% Solution z3(t) returned on the time interval [t03 tf3]

v3rd	   	= z3(:,1);	
A3rd  		= z3(:,2);	
phi3rd	   	= z3(:,3);						
r3rd	   	= z3(:,4);						
delta3rd	= z3(:,5);						
lamda3rd	= z3(:,6);
m3rd		= z3(:,7);

% Initial conditions for the fourth flight phase (ballistic)

v4th_0		= v3rd(end);					
A4th_0		= A3rd(end);				
phi4th_0	= phi3rd(end);					
r4th_0		= r3rd(end);					
delta4th_0	= delta3rd(end);					
lamda4th_0	= lamda3rd(end);
m4th_0		= m0_nikha+m_pay;

t04	= tf3;
tf4	= t04+2; % Two seconds between the stages
tspan4	= [t04 tf4];

z4_0		= [v4th_0 A4th_0 phi4th_0 r4th_0 delta4th_0 lamda4th_0 m4th_0];

[t4,z4]=ode45(@ballistic,tspan4,z4_0,options1);

v4th	   	= z4(:,1);						
A4th  		= z4(:,2);						
phi4th	   	= z4(:,3);						
r4th	   	= z4(:,4);						
delta4th	= z4(:,5);						
lamda4th	= z4(:,6);
m4th		= z4(:,7);

% ...Initial conditions for the fifth flight phase (powered)

v5th_0		= v4th(end);					
A5th_0		= A4th(end);				
phi5th_0	= phi4th(end);					
r5th_0		= r4th(end);					
delta5th_0	= delta4th(end);				
lamda5th_0	= lamda4th(end);
m5th_0		= m4th(end);

t05	= tf4;
tf5	= t05+tburn_nikha;
tspan5	= [t05 tf5];

z5_0		= [v5th_0 A5th_0 phi5th_0 r5th_0 delta5th_0 lamda5th_0 m5th_0];

[t5,z5]=ode45(@powered_nikha,tspan5,z5_0,options1);

v5th	   	= z5(:,1);						
A5th  		= z5(:,2);						
phi5th	   	= z5(:,3);						
r5th	   	= z5(:,4);						
delta5th	= z5(:,5);						
lamda5th	= z5(:,6);
m5th		= z5(:,7);

% ...Initial conditions for the fifth flight phase (3rd powered)

v6th_0		= v5th(end);					
A6th_0		= A5th(end);				
phi6th_0	= phi5th(end);					
r6th_0		= r5th(end);					
delta6th_0	= delta5th(end);					
lamda6th_0	= lamda5th(end);
m6th_0		= m_pay;

t06	= tf5;
tf6	= t06+10000;   % Enough time to allow the capsule reachs the apogee point 
tspan6	= [t06 tf6];

z6_0		= [v6th_0 A6th_0 phi6th_0 r6th_0 delta6th_0 lamda6th_0 m6th_0];

[t6,z6]=ode45(@ballistic,tspan6,z6_0,options1);

v6th	   	= z6(:,1);						
A6th  		= z6(:,2);						
phi6th	   	= z6(:,3);						
r6th	   	= z6(:,4);						
delta6th	= z6(:,5);						
lamda6th	= z6(:,6);
m6th		= z6(:,7);

% Initial conditions

v7th_0		= v6th(end);
A7th_0		= A6th(end);
phi7th_0	= phi6th(end);
r7th_0		= r6th(end);
delta7th_0	= delta6th(end);
lamda7th_0	= lamda6th(end);
m7th_0		= m6th(end);

t07		= t6(end);
tf7		= t07+10000;
tspan7		= [t07 tf7];

z7_0		= [v7th_0 A7th_0 phi7th_0 r7th_0 delta7th_0 lamda7th_0 m7th_0];

[t7,z7]=ode45(@reentry,tspan7,z7_0,options2);

v7th	   	= z7(:,1);						
A7th  		= z7(:,2);						
phi7th	   	= z7(:,3);						
r7th	   	= z7(:,4);						
delta7th	= z7(:,5);						
lamda7th	= z7(:,6);
m7th		= z7(:,7);


% ...Total flight data vectors

t 	 = 		[t1; t2; t3; t4; t5; t6; t7];
v 	 = 		[v1st; v2nd; v3rd; v4th; v5th; v6th; v7th]*1.e-3;
A	 = 		[A1st; A2nd; A3rd; A4th; A5th; A6th; A7th]/deg;
phi 	 =		[phi1st; phi2nd; phi3rd; phi4th; phi5th; phi6th; phi7th]/deg;
r	 =		[r1st; r2nd; r3rd; r4th; r5th; r6th; r7th]*1.e-3;
delta 	 = 		[delta1st; delta2nd; delta3rd; delta4th; delta5th; delta6th; delta7th]/deg;
del	 =		[delta1st; delta2nd; delta3rd; delta4th; delta5th; delta6th; delta7th];
lamda	 =		[lamda1st; lamda2nd; lamda3rd; lamda4th; lamda5th; lamda6th; lamda7th]/deg;
m 	 = 		[m1st; m2nd; m3rd; m4th; m5th; m6th; m7th];

R	 =		Re*(1-epsilon*(sin(del)).^2)*1.e-3;

h	 =		abs(r-R);

h_100 = find(h>100);	% Flight in space, over 100 kilometers

over100_0 = h_100(1);
over100_f = h_100(end);

microg_time = t(over100_f) - t(over100_0);

Q = dynamic_pressure(v,h); % Call the dynamic pressure function which gives the dynamic pressure [kPa]

figure(1)
plot(t,h)
axis equal
xlabel('Time (s)')
ylabel('Altitude (km)')
axis([-inf, inf, -inf, inf])
grid

figure(2)
plot(v,h,'m')
xlabel('Speed (km/s)')
ylabel('Altitude (km)')
axis([-inf, inf, -inf, inf])
grid

figure(3)
plot(t, v)
xlabel('time (s)')
ylabel('Speed (km/s)')
grid

figure(4)
plot(t,phi)
xlabel ('Time')
ylabel ('Angle (deg)')
grid

figure(5)
plot(t,m)
xlabel('Time (s)')
ylabel('Mass (Kg)')
axis([-inf, 100, -inf, inf])
grid

figure(6)
plot(lamda,delta,'r')
xlabel('Longitude')
ylabel('Latitude')
grid

figure(7)
plot(h,v)
xlabel('Altitude (km)')
ylabel('velocity (km/s)')
grid

figure(8)
plot(Q,h,'r')
xlabel('Dynamic pressure (kPa)')
ylabel('Altitude (Km)')
axis([-inf, inf, -inf, 80])
grid



