function xprime = ballistic(t,x);

global gamma g0 rho0 hscale Re Ma15C subs trans a b c d e f A_talos Th_talos mdot_talos omega epsilon Lc_fustrum A_fustrum

v=x(1);
A=x(2);
phi=x(3);
r=x(4);
delta=x(5);
lamda=x(6);
m=x(7);

R = Re*(1-epsilon*(sin(delta))^2);
h = r-R;

Lc = Lc_fustrum;
Y = atmosphere(h,Lc);
rho = Y(1);
v_sound = Y(2);
Knudsen = Y(3);
Ma = v/v_sound;
s = Ma*sqrt(gamma/2);
continuum_flow = 0.01;
molecular_flow = 10;

Cd = 0.4;

[g_center g_delta]=gravity(r,delta);

gc=abs(-g_center); % This action reverses the gravity vector sense.
gdelta=-g_delta;


D = 1/2*rho*v^2*A_fustrum*Cd;

xprime = [-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-omega^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)); ...
0; ...
0; ...
v*sin(phi); ...
v/r*cos(phi)*cos(A); ...
v*cos(phi)*sin(A)/r*cos(delta); ...
0; ...
-D/m; ...
-abs(-gc*sin(phi)+gdelta*cos(phi)*cos(A))/m];
