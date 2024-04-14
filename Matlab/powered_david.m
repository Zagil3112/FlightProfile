function xprime = powered(t,x,p)

global gamma Re subs trans a b c d e f omega epsilon Lc_fustrum A_fustrum



v=x(1);
A=x(2);
phi=x(3);
r=x(4);
delta=x(5);
lamda=x(6);
m=x(7);

A_rocket=p(1);
Th_rocket=p(2);
mdot_rocket=p(3);

R = Re*(1-epsilon*(sin(delta))^2);

Lc = Lc_fustrum;
h = r-R;
Y = atmosphere(h,Lc);
rho = Y(1);
v_sound = Y(2);
Knudsen = Y(3);
Ma = v/v_sound;
s = Ma*sqrt(gamma/2);
continuum_flow = 0.01;
molecular_flow = 10;

Cd=0.4;
   


[g_center g_delta]=gravity(r,delta);

gc=abs(-g_center); % This action reverses the gravity vector sense.
gdelta=-g_delta;

D = 1/2*rho*v^2*A_rocket*Cd;


xprime = [Th_rocket/m-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-omega^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)); ...
0; ...
0; ...
v*sin(phi); ...
v/r*cos(phi)*cos(A); ...
v*cos(phi)*sin(A)/r*cos(delta); ...
mdot_rocket; ...
-D/m; ...
-abs(-gc*sin(phi)/m)];
