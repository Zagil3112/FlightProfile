% ...THIS PROGRAM CALCULATES THE TRAJECTORY FOR A THREE STAGE SOUNDING ROCKET

% DECLARE GLOBAL VARIABLES FOR USE IN EVERY FUNCTION
global rho0 gamma R_air Re R_polar subs trans a b c d e f omega diam_talos Lc_capsule A_capsule Lc_fustrum A_fustrum epsilon mu J2 J3 J4

deg	= pi/180;			% Convert degrees to radians
rho0	= 1.225;				% Sea level atmosphere density
gamma	= 1.4;				% Sea level specific heat ratio
R_air	= 257.058;			% Sea level gas constant for air [J/Kg-K]

% NGA/NASA EGM96, N=M=360, WGS84 REVISED EARTH GRAVITATIONAL MODEL
Re	= 6378.137e3;						% Semi-major axis of WGS 84 Ellipsoid (Radius of the earth equator for the NGA/NASA EGM96 WGS84) Gravitational model[m]
R_polar	= 6356.752314245e3;		% Semi-minor axis of WGS 84 Ellipsoid (Polar radius of the earth) [m]
epsilon	= 1.0/298.2572235630;	% Flattening of WGS 84 Ellipsoid
omega	= 7.2921150e-5;				% Angular speed of Earth's rotation [rad/sec]
mu	= 3.986004418e14; 				% mu=GMe
J2	= 1.08263e-3;						% J2, J3, J4 harmonics for modelling the WGS84
J3	= 2.532153e-7;
J4	= 1.6109876e-7;
g = 9.80665; 							% standard Earth gravity

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
diam_talos	= 0.76;									% Vehicle Diameter [m]
A_talos		= pi/4*(diam_talos)^2;				% Frontal area [m^2]
m0_talos	= 1996;										% Talos total mass [kg]
mfinal_talos	= 496;								% Talos Unfueled mass [kg]
Th_talos 	= 516e3;									% Talos thrust [N]
mdot_talos	=-288.46;								% Talos propellant consumption rate [kg/s]
mprop_talos	= m0_talos - mfinal_talos;			% Talos Propellant mass [kg]
tburn_talos	= mprop_talos/abs(mdot_talos);	% Talos burning time [s]

% ROCKET CONSTANTS FOR TOMAHAWK STAGE
diam_toma	= 0.58;									% Vehicle Diameter [m]
A_toma		= pi/4*(diam_toma)^2;				% Frontal area [m^2]
m0_toma		= 1363;									% Tomahawk total mass [kg]
mfinal_toma	= 602;									% Tomahawk unfueled mass
Th_toma		= 457e3;									% Tomahawk thrust [N]
mdot_toma	= -217.43;								% Tomahawk propellant consumption rate [kg/s]
mprop_toma	= m0_toma - mfinal_toma;			% Tomahawk propellant mass [kg]
tburn_toma	= mprop_toma/abs(mdot_toma);		% Tomahawk burning time [s]

% ROCKET CONSTANTS FOR NIKHA STAGE
diam_nikha		= 0.44;									% Vehicle Diameter [m]
A_nikha			= pi/4*(diam_nikha)^2;				% Frontal area [m^2]
m0_nikha			= 399;									% Nikha total mass [kg]
mfinal_nikha	= 70;										% Nikha unfueled mass [kg]
Th_nikha			= 50500;									% Nikha thrust [N]
mdot_nikha		= -19.35;								% Nikha propellant consumption rate [Kg/s]
mprop_nikha		= m0_nikha - mfinal_nikha;			% Nikha propellant mass [kg]
tburn_nikha		= mprop_nikha/abs(mdot_nikha);	% Nikha burning time [s]

% CONSTANTS FOR ROCKET FUSTRUM
diam_fustrum	= 0.8;							% Fustrum diameter [m]
A_fustrum		= pi/4*(diam_fustrum)^2;	% Fustrum area [m^2]
Lc_fustrum		= 0.2;							% Characteristic longitude for ascent Knudsen [m]

% CONSTANTS FOR CAPSULE
diam_capsule	= 0.78;							% Reentry capsule diameter [m]
A_capsule	= pi/4*(diam_capsule)^2;		% Reentry capsule frontal area [m]
Lc_capsule	= 0.21;								% Characteristic longitude for reentry Knudsen [m]

% PAYLOAD MASS
m_pay		= 120; % [kg]

% INITIAL CONDITIONS
t0 = 0;					% The flight time clock begins
tf1 = tburn_talos;	% Talos burning time for time interval
tspan1 = [t0 tf1];	% Time interval for first flight phase integration 

v1st_0		= 0;												% Initial velocity
A1st_0		= 90*deg;										% Initial azimuth
phi1st_0	= 85*deg;											% Initial angle wrt local horizontal
delta1st_0	= -37.940194*deg;								% Launch site latitude
lamda1st_0	= -75.466389*deg;								% Launch site longitude
r1st_0		= Re*(1-epsilon*(sin(delta1st_0))^2);	% Initial datum sea level due to WGS84 model
m1st_0		= m0_talos+m0_toma+m0_nikha+m_pay;		% Total initial mass
vlossdrag1st_0	= 0;											% Initial velocity loss	due to drag
vlossgrav1st_0	= 0;											% Initial velocity loss due to gravity

z1_0=[v1st_0 A1st_0 phi1st_0 r1st_0 delta1st_0 lamda1st_0 m1st_0 vlossdrag1st_0 vlossgrav1st_0];

% Specify tolerance and event stop at apogee for ode45 integrator
options1 = odeset('Events',@ballis_event,'reltol',1e-12,'abstol',1e-12);
% Specify tolerance and event stop at the end of the flight
options2 = odeset('Events',@reentry_event,'reltol',1e-12,'abstol',1e-12);

p1=[A_talos Th_talos mdot_talos]; % Talos rocket parameters to be sent to powered.m

[t1,z1]=ode45(@powered,tspan1,z1_0,[],p1);

% Solution z(t) returned on the time interval [t0 tf] for first flight phase:
v1st				= z1(:,1);	 % Velocity [m/s]
A1st				= z1(:,2);	 % Azimuth angle [rad]
phi1st			= z1(:,3);	 % Flight path angle [rad]
r1st				= z1(:,4);	 % Distance from Earth center [m]
delta1st			= z1(:,5);	 % Latitude [rad]
lamda1st			= z1(:,6);	 % Longitude [rad]
m1st				= z1(:,7);	 % Mass [kg]
vlossdrag1st	= z1(:,8);	 % Velocity loss due to drag [m/s]
vlossgrav1st	= z1(:,9);	 % Velocity loss due to gravity [m/s]

% Initial conditions for the second flight phase (1st ballistic)
% the same quantities as in  flight phase 1 apply

v2nd_0			= v1st(end);				
A2nd_0			= A1st(end);				
phi2nd_0			= phi1st(end);					
r2nd_0			= r1st(end);					
delta2nd_0		= delta1st(end);					
lamda2nd_0		= lamda1st(end);
m2nd_0			= m0_toma+m0_nikha+m_pay; 	% Talos stage has been detached
vlossdrag2nd_0	= vlossdrag1st(end);
vlossgrav2nd_0	= vlossgrav1st(end);

t02	 = tf1;			% Initial time for 2nd phase = final time 1st phase
tf2	 = t02+2;		% Two seconds between the stages	
tspan2 = [t02, tf2];	% Time interval for second flight phase integration

z2_0 		= [v2nd_0 A2nd_0 phi2nd_0 r2nd_0 delta2nd_0 lamda2nd_0 m2nd_0 vlossdrag2nd_0 vlossgrav2nd_0];

[t2,z2]=ode45(@ballistic,tspan2,z2_0);

% Solution z2(t) returned on the time interval [t02 tf2]
v2nd	   	 = z2(:,1);	% Velocity [m/s]
A2nd  		 = z2(:,2);	% Azimuth angle [rad]
phi2nd	    = z2(:,3);	% Flight path angle [rad]
r2nd	   	 = z2(:,4);	% Distance from Earth Center [m]
delta2nd		 = z2(:,5);	% Latitude [rad]
lamda2nd		 = z2(:,6);	% Longitude [rad]
m2nd			 = z2(:,7);	% Mass [kg]
vlossdrag2nd = z2(:,8);	% Velocity loss due to drag [m/s]
vlossgrav2nd = z2(:,9);	% Velocity loss due to gravity [m/s]

% Initial conditions for the third flight phase (powered)
% % the same quantities as in  flight phase 1 and 2 apply

v3rd_0			= v2nd(end);					
A3rd_0			= A2nd(end);				
phi3rd_0			= phi2nd(end);					
r3rd_0			= r2nd(end);					
delta3rd_0		= delta2nd(end);					
lamda3rd_0		= lamda2nd(end);
m3rd_0			= m2nd(end);
vlossdrag3rd_0	= vlossdrag2nd(end);
vlossgrav3rd_0	= vlossgrav2nd(end);

t03		= tf2;  % Initial time for 3rd phase = 2nd phase final time
tf3		= t03+tburn_toma;
tspan3	= [t03 tf3];

z3_0		= [v3rd_0 A3rd_0 phi3rd_0 r3rd_0 delta3rd_0 lamda3rd_0 m3rd_0 vlossdrag3rd_0 vlossgrav3rd_0];

p2=[A_toma Th_toma mdot_toma]; % Tomahawk rocket parameters to be sent to powered.m

[t3,z3]=ode45(@powered,tspan3,z3_0,[],p2);

% Solution z3(t) returned on the time interval [t03 tf3]
% The same quantities as for previous phases apply

v3rd	   		= z3(:,1);	
A3rd  			= z3(:,2);	
phi3rd	   	= z3(:,3);						
r3rd	   		= z3(:,4);						
delta3rd			= z3(:,5);						
lamda3rd			= z3(:,6);
m3rd				= z3(:,7);
vlossdrag3rd	= z3(:,8);
vlossgrav3rd	= z3(:,9);

% Initial conditions for the fourth flight phase (ballistic)
% The same quantities as for previous phases apply

v4th_0			= v3rd(end);					
A4th_0			= A3rd(end);				
phi4th_0			= phi3rd(end);					
r4th_0			= r3rd(end);					
delta4th_0		= delta3rd(end);					
lamda4th_0		= lamda3rd(end);
m4th_0			= m0_nikha+m_pay;
vlossdrag4th_0	= vlossdrag3rd(end);
vlossgrav4th_0	= vlossgrav3rd(end);

t04	= tf3;	% Initial time for 4th phase = 3rd phase final time
tf4	= t04+2; % Two seconds between the stages
tspan4	= [t04 tf4];

z4_0		= [v4th_0 A4th_0 phi4th_0 r4th_0 delta4th_0 lamda4th_0 m4th_0 vlossdrag4th_0 vlossgrav4th_0];

[t4,z4]=ode45(@ballistic,tspan4,z4_0);

v4th	   	 = z4(:,1);						
A4th  		 = z4(:,2);						
phi4th	    = z4(:,3);						
r4th	   	 = z4(:,4);						
delta4th		 = z4(:,5);						
lamda4th		 = z4(:,6);
m4th			 = z4(:,7);
vlossdrag4th = z4(:,8);
vlossgrav4th = z4(:,9);

% Initial conditions for the fifth flight phase (Nikha-powered)
% The same quantities as for previous phases apply

v5th_0		   = v4th(end);					
A5th_0		   = A4th(end);				
phi5th_0		   = phi4th(end);					
r5th_0		   = r4th(end);					
delta5th_0	   = delta4th(end);				
lamda5th_0	   = lamda4th(end);
m5th_0		   = m4th(end);
vlossdrag5th_0	= vlossdrag4th(end);
vlossgrav5th_0	= vlossgrav4th(end);

t05	 = tf4;
tf5	 = t05+tburn_nikha;
tspan5 = [t05 tf5];

z5_0	 = [v5th_0 A5th_0 phi5th_0 r5th_0 delta5th_0 lamda5th_0 m5th_0 vlossdrag5th_0 vlossgrav5th_0];

p3=[A_nikha Th_nikha mdot_nikha]; % Nikha rocket parameters to be sent to powered.m function

[t5,z5]=ode45(@powered,tspan5,z5_0,[],p3);

v5th	   	 = z5(:,1);						
A5th  		 = z5(:,2);						
phi5th	    = z5(:,3);						
r5th	   	 = z5(:,4);						
delta5th		 = z5(:,5);						
lamda5th		 = z5(:,6);
m5th			 = z5(:,7);
vlossdrag5th = z5(:,8);
vlossgrav5th = z5(:,9);

% Initial conditions for the sixth flight phase (3rd ballistic)
% The same quantities as for previous phases apply

v6th_0			= v5th(end);					
A6th_0			= A5th(end);				
phi6th_0			= phi5th(end);					
r6th_0			= r5th(end);					
delta6th_0		= delta5th(end);					
lamda6th_0		= lamda5th(end);
m6th_0			= m_pay;  % The nikha stage and the fustrum has been detached, only the capsule keeps on flight
vlossdrag6th_0	= vlossdrag5th(end);
vlossgrav6th_0	= vlossgrav5th(end);

t06	= tf5;
tf6	= t06+10000;    % Enough time to allow the capsule reachs the apogee point. 
							 % The integration will be stopped when ode45 reachs the event
							 % stop at apogee written at ballis_event.m
tspan6	= [t06 tf6];

z6_0		= [v6th_0 A6th_0 phi6th_0 r6th_0 delta6th_0 lamda6th_0 m6th_0 vlossdrag6th_0 vlossgrav6th_0];

[t6,z6]=ode45(@ballistic,tspan6,z6_0,options1);

v6th	       = z6(:,1);						
A6th  	    = z6(:,2);						
phi6th	    = z6(:,3);						
r6th	       = z6(:,4);						
delta6th	    = z6(:,5);						
lamda6th	    = z6(:,6);
m6th		    = z6(:,7);
vlossdrag6th = z6(:,8);
vlossgrav6th = z6(:,9);

% Initial conditions for the seventh flight phase (reentry)
% The same quantities as for previous phases apply

v7th_0		= v6th(end);
A7th_0		= A6th(end);
phi7th_0		= phi6th(end);
r7th_0		= r6th(end);
delta7th_0	= delta6th(end);
lamda7th_0	= lamda6th(end);
m7th_0		= m6th(end);

t07	 = t6(end);
tf7	 = t07+10000;   % Enough time to allow the capsule reachs the final flight point. 
							 % The integration will be stopped when ode45 reachs the event
							 % stop at 5 km. altitude over sea level, written at reentry_event.m
tspan7 = [t07 tf7];

z7_0 = [v7th_0 A7th_0 phi7th_0 r7th_0 delta7th_0 lamda7th_0 m7th_0];

[t7,z7]=ode45(@reentry,tspan7,z7_0,options2);

v7th	   = z7(:,1);						
A7th  	= z7(:,2);						
phi7th	= z7(:,3);						
r7th	   = z7(:,4);						
delta7th	= z7(:,5);						
lamda7th	= z7(:,6);
m7th		= z7(:,7);

% TOTAL FLIGHT DATA VECTORS

t 	 		 = 	[t1; t2; t3; t4; t5; t6; t7];
v 	 		 = 	[v1st; v2nd; v3rd; v4th; v5th; v6th; v7th]*1.e-3;
v_asc		 =    [v1st; v2nd; v3rd; v4th; v5th; v6th]*1.e-3;
v_des		 =		[v7th]*1.e-3;
A	 		 = 	[A1st; A2nd; A3rd; A4th; A5th; A6th; A7th]/deg;
phi 	 	 =		[phi1st; phi2nd; phi3rd; phi4th; phi5th; phi6th; phi7th]/deg;
r	 		 =		[r1st; r2nd; r3rd; r4th; r5th; r6th; r7th]*1.e-3;
delta 	 = 	[delta1st; delta2nd; delta3rd; delta4th; delta5th; delta6th; delta7th]/deg;
del	 	 =		[delta1st; delta2nd; delta3rd; delta4th; delta5th; delta6th; delta7th];
lamda	 	 =		[lamda1st; lamda2nd; lamda3rd; lamda4th; lamda5th; lamda6th; lamda7th]/deg;
m 	 		 = 	[m1st; m2nd; m3rd; m4th; m5th; m6th; m7th];
vlossdrag =		[vlossdrag1st; vlossdrag2nd; vlossdrag3rd; vlossdrag4th; vlossdrag5th; vlossdrag6th]*1.e-3;
vlossgrav =		[vlossgrav1st; vlossgrav2nd; vlossgrav3rd; vlossgrav4th; vlossgrav5th; vlossgrav6th]*1.e-3;

R1st	 	 =		Re*(1-epsilon*(sin(delta1st)).^2);
R2nd	 	 =		Re*(1-epsilon*(sin(delta2nd)).^2);
R3rd	 	 =		Re*(1-epsilon*(sin(delta3rd)).^2);
R4th	 	 =		Re*(1-epsilon*(sin(delta4th)).^2);
R5th	 	 =		Re*(1-epsilon*(sin(delta5th)).^2);
R6th	 	 =		Re*(1-epsilon*(sin(delta6th)).^2);
R7th	 	 =		Re*(1-epsilon*(sin(delta7th)).^2);

h1st	 	 =		abs((r1st-R1st)*1.e-3);
h2nd	 	 =		abs((r2nd-R2nd)*1.e-3);
h3rd	 	 =		abs((r3rd-R3rd)*1.e-3);
h4th	 	 =		abs((r4th-R4th)*1.e-3);
h5th	 	 =		abs((r5th-R5th)*1.e-3);
h6th	 	 =		abs((r6th-R6th)*1.e-3);
h7th	 	 =		abs((r7th-R7th)*1.e-3);
h	 		 =		[h1st; h2nd; h3rd; h4th; h5th; h6th; h7th];
h_asc		 =    [h1st; h2nd; h3rd; h4th; h5th; h6th];

% RANGE ESTIMATION FROM VINCENTY FORMULA

Ran1st	 =		distVincenty(delta1st(1),delta1st(end),lamda1st(1),lamda1st(end));
Ran2nd	 =		distVincenty(delta1st(1),delta2nd(end),lamda1st(1),lamda2nd(end));
Ran3rd	 =		distVincenty(delta1st(1),delta3rd(end),lamda1st(1),lamda3rd(end));
Ran4th	 =		distVincenty(delta1st(1),delta4th(end),lamda1st(1),lamda4th(end));
Ran5th	 =		distVincenty(delta1st(1),delta5th(end),lamda1st(1),lamda5th(end));
Ran6th	 =		distVincenty(delta1st(1),delta6th(end),lamda1st(1),lamda6th(end));
Ran7th	 =		distVincenty(delta1st(1),delta7th(end),lamda1st(1),lamda7th(end));
RanApo	 =		distVincenty(delta1st(1),delta7th(1),lamda1st(1),lamda7th(1));

% G ACCELERATIONS
%flight phase 1=ascent; 2=reentry;
ascent	 = 	1; 
reent		 =		2;
coast		 =		0; % if there's not thrust
%constants
p1st		 =		[Th_talos A_fustrum Lc_fustrum ascent];
p2nd		 =		[coast A_fustrum Lc_fustrum ascent];
p3rd		 = 	[Th_toma A_fustrum Lc_fustrum ascent];
p4th		 =		[coast A_fustrum Lc_fustrum ascent];
p5th		 =		[Th_nikha A_fustrum Lc_fustrum ascent];
p6th		 = 	[coast A_fustrum Lc_fustrum ascent];
p7th		 = 	[coast A_capsule Lc_capsule reent];

Gs1st		 =		gacceler(p1st,v1st,A1st,phi1st,r1st,delta1st,m1st);
Gs2nd		 =		gacceler(p2nd,v2nd,A2nd,phi2nd,r2nd,delta2nd,m2nd);
Gs3rd		 =		gacceler(p3rd,v3rd,A3rd,phi3rd,r3rd,delta3rd,m3rd);
Gs4th	    =		gacceler(p4th,v4th,A4th,phi4th,r4th,delta4th,m4th);
Gs5th	    =		gacceler(p5th,v5th,A5th,phi5th,r5th,delta5th,m5th);
Gs6th	    =		gacceler(p6th,v6th,A6th,phi6th,r6th,delta6th,m6th);
Gs7th	    =		gacceler(p7th,v7th,A7th,phi7th,r7th,delta7th,m7th);

Gs1st_xv=Gs1st(:,1); Gs1st_yv=Gs1st(:,2); Gs1st_zv=Gs1st(:,3); Gs1st_total=Gs1st(:,4);
Gs2nd_xv=Gs2nd(:,1); Gs2nd_yv=Gs2nd(:,2); Gs2nd_zv=Gs2nd(:,3); Gs2nd_total=Gs2nd(:,4);
Gs3rd_xv=Gs3rd(:,1); Gs3rd_yv=Gs3rd(:,2); Gs3rd_zv=Gs3rd(:,3); Gs3rd_total=Gs3rd(:,4);
Gs4th_xv=Gs4th(:,1); Gs4th_yv=Gs4th(:,2); Gs4th_zv=Gs4th(:,3); Gs4th_total=Gs4th(:,4);
Gs5th_xv=Gs5th(:,1); Gs5th_yv=Gs5th(:,2); Gs5th_zv=Gs5th(:,3); Gs5th_total=Gs5th(:,4);
Gs6th_xv=Gs6th(:,1); Gs6th_yv=Gs6th(:,2); Gs6th_zv=Gs6th(:,3); Gs6th_total=Gs6th(:,4);
Gs7th_xv=Gs7th(:,1); Gs7th_yv=Gs7th(:,2); Gs7th_zv=Gs7th(:,3); Gs7th_total=Gs7th(:,4);

Gxv=[Gs1st_xv; Gs2nd_xv; Gs3rd_xv; Gs4th_xv; Gs5th_xv; Gs6th_xv; Gs7th_xv]/g;
Gyv=[Gs1st_yv; Gs2nd_yv; Gs3rd_yv; Gs4th_yv; Gs5th_yv; Gs6th_yv; Gs7th_yv]/g;
Gzv=[Gs1st_zv; Gs2nd_zv; Gs3rd_zv; Gs4th_zv; Gs5th_zv; Gs6th_zv; Gs7th_zv]/g;
Gtotal=[Gs1st_total; Gs2nd_total; Gs3rd_total; Gs4th_total; Gs5th_total; Gs6th_total; Gs7th_total]/g;


h_100 = find(h>100);	% Flight in space, over 100 kilometers

over100_0 = h_100(1);
over100_f = h_100(end);

microg_time = t(over100_f) - t(over100_0);

Q = dynamic_pressure(v,h); % Call the dynamic pressure function which gives the dynamic pressure [kPa]
Q_asc = dynamic_pressure(v_asc,h_asc);
Q_des = dynamic_pressure(v_des,h7th);

fprintf('\n --------------FIRST FLIGHT PHASE-----------------------')
fprintf('\n                     Powered                            ')
fprintf('\n Final speed at Talos burnout        = %10g Km/s',v1st(end)*1.e-3)
fprintf('\n Altitude at Talos burnout           = %10g Km',h1st(end))
fprintf('\n Range at Talos burnout              = %10g Km',Ran1st*1.e-3)
fprintf('\n Angle at Talos burnout              = %10g deg',phi1st(end)/deg)
fprintf('\n Mass at Talos burnout               = %10g Kg',m1st(end))
fprintf('\n Propellant mass burned              = %10g Kg',mprop_talos)
fprintf('\n Talos ignition time                 = %10g s',tburn_talos)
fprintf('\n Mission elapsed time                = %10g s',t1(end))
fprintf('\n')
fprintf('\n --------------SECOND FLIGHT PHASE----------------------')
fprintf('\n                      Ballistic                         ')
fprintf('\n Mass at 1st Ballistic flight phase  = %10g Kg',m2nd(1))
fprintf('\n Final speed at 1st Ballistic phase  = %10g Km/s',v2nd(end)*1.e-3)
fprintf('\n Final altitude at 1st Ballis phase  = %10g Km',h2nd(end))
fprintf('\n Range after 1st Ballistic phase     = %10g Km',Ran2nd*1.e-3)
fprintf('\n Angle after 1st Ballistic phase     = %10g deg',phi2nd(end)/deg)
fprintf('\n 1st Ballistic flight time           = %10g s',t2(end)-t2(1))
fprintf('\n Mission elapsed time                = %10g s',t2(end))
fprintf('\n')
fprintf('\n --------------THIRD FLIGHT PHASE-----------------------')
fprintf('\n                     Powered                            ')
fprintf('\n Final speed at Tomahawk burnout     = %10g Km/s',v3rd(end)*1.e-3)
fprintf('\n Altitude at Tomahawk burnout        = %10g Km',h3rd(end))
fprintf('\n Range at Tomahawk burnout           = %10g Km',Ran3rd*1.e-3)
fprintf('\n Angle at Tomahawk burnout           = %10g deg',phi3rd(end)/deg)
fprintf('\n Mass at Tomahawk burnout            = %10g Kg',m3rd(end))
fprintf('\n Propellant mass burned              = %10g Kg',mprop_talos)
fprintf('\n Tomahawk ignition time              = %10g s',tburn_toma)
fprintf('\n Mission elapsed time                = %10g s',t3(end))
fprintf('\n')
fprintf('\n --------------FOURTH FLIGHT PHASE----------------------')
fprintf('\n                     Ballistic                          ')
fprintf('\n Mass at 2nd Ballistic flight phase  = %10g Kg',m4th(1))
fprintf('\n Final speed at 2nd Ballistic phase  = %10g Km/s',v4th(end)*1.e-3)
fprintf('\n Final altitude at 2nd Ballis phase  = %10g Km',h4th(end))
fprintf('\n Range after 2nd Ballistic phase     = %10g Km',Ran4th*1.e-3)
fprintf('\n Angle after 2nd Ballistic phase     = %10g deg',phi4th(end)/deg)
fprintf('\n 2nd Ballistic flight time           = %10g s',t4(end)-t4(1))
fprintf('\n Mission elapsed time                = %10g s',t4(end))
fprintf('\n')
fprintf('\n --------------FIFTH FLIGHT PHASE-----------------------')
fprintf('\n                     Powered                            ')
fprintf('\n Final speed at Nikha burnout        = %10g Km/s',v5th(end)*1.e-3)
fprintf('\n Altitude at Nikha burnout           = %10g Km',h5th(end))
fprintf('\n Range at Nikha burnout              = %10g Km',Ran5th*1.e-3)
fprintf('\n Angle at Nikha burnout              = %10g deg',phi5th(end)/deg)
fprintf('\n Mass at Nikha burnout               = %10g Kg',m5th(end))
fprintf('\n Propellant mass burned              = %10g Kg',mprop_nikha)
fprintf('\n Nikha ignition time                 = %10g s',tburn_nikha)
fprintf('\n Mission elapsed time                = %10g s',t5(end))
fprintf('\n')
fprintf('\n --------------SIXTH FLIGHT PHASE-----------------------')
fprintf('\n                     Ballistic                          ')
fprintf('\n Mass at 3rd Ballistic flight phase  = %10g Kg',m6th(1))
fprintf('\n Final speed at 3rd Ballistic phase  = %10g Km/s',v6th(end)*1.e-3)
fprintf('\n Final altitude at 3rd Ballis phase  = %10g Km',h6th(end))
fprintf('\n Range after 3rd Ballistic phase     = %10g Km',Ran6th*1.e-3)
fprintf('\n Angle after 3rd Ballistic phase     = %10g deg',phi6th(end)/deg)
fprintf('\n 3rd Ballistic flight time           = %10g s',t6(end)-t6(1))
fprintf('\n Mission elapsed time                = %10g s',t6(end))
fprintf('\n')
fprintf('\n ---------------SEVENTH FLIGHT PHASE--------------------')
fprintf('\n                     Reentry                            ')
fprintf('\n Final speed                         = %10g Km/s',v(end))
fprintf('\n Speed at apogee                     = %10g Km/s',v7th(1)*1.e-3)
fprintf('\n Apogee Altitude                     = %10g Km',max(h))
fprintf('\n Apogee Range                        = %10g Km',RanApo*1.e-3)
fprintf('\n Total Range                         = %10g Km',Ran7th*1.e-3)
fprintf('\n Final angle                         = %10g deg',phi(end))
fprintf('\n Final Mass                          = %10g Kg',m(end))
fprintf('\n Reentry flight time                 = %10g s',t7(end)-t7(1))
fprintf('\n Microgravity time                   = %10g s',microg_time)
fprintf('\n Total flight time                   = %10g s',t(end))
fprintf('\n Final altitude                      = %10g Km',h(end))
fprintf('\n Drag loss (only ascent)             = %10g Km/s',vlossdrag(end))
fprintf('\n Gravity loss (only ascent)          = %10g Km/s',vlossgrav(end))
fprintf('\n\n ---------------------------------------------------\n')

figure(1)
plot(t,h)
axis equal
xlabel('Time (s)')
ylabel('Altitude (km)')
axis([-inf, inf, -inf, inf])
grid

figure(2)
plot(t, v)
xlabel('time (s)')
ylabel('velocity (km/s)')
grid

figure(3)
plot(t,phi)
xlabel ('Time')
ylabel ('Angle (deg)')
grid

figure(4)
plot(t,m)
xlabel('Time (s)')
ylabel('Mass (Kg)')
axis([0, 100, 0, inf])
grid

figure(5)
plot(lamda,delta,'g')
xlabel('Longitude')
ylabel('Latitude')

figure(6)
plot(h_asc,v_asc,'r',h7th,v_des,'b')
xlabel('altitude (km)')
ylabel('velocity (km/s)')
grid

figure(7)
plot(Q_asc,h_asc,'r',Q_des,h7th,'g')
xlabel('Dynamic pressure (kPa)')
ylabel('Altitude (Km)')
axis([-inf, inf, -inf, 80])
grid

figure(8)
plot(t,Q,'g')
xlabel('Time')
ylabel('Dynamic pressure ascent (kPa)')
axis([-inf, 60, -inf, inf])
grid

figure(9)
plot(t,Gtotal,'r')
xlabel('Time')
ylabel('Ascent Gs')
axis([-inf, 35, -inf, inf])
grid



