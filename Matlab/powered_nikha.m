function xprime = powered_nikha(t,x);

global g0 rho0 hscale Re Ma15C subs trans a b c d e f A_nikha Th_nikha mdot_nikha omega epsilon

v=x(1);
A=x(2);
phi=x(3);
r=x(4);
delta=x(5);
lamda=x(6);
m=x(7);

R = Re*(1-epsilon*(sin(delta))^2);
h = r-R;
Y = atmosphere(h);
rho = Y(1);
v_sound = Y(2);
Ma = v/v_sound;

if Ma < subs
   Cd = a;
   elseif Ma < trans
	Cd = -b+c*Ma;
    else
	Cd = d+e/sqrt(Ma^2-f);
end

[g_center g_delta]=gravity(r,delta);

gc=abs(-g_center); % This action reverses the gravity vector sense.
gdelta=-g_delta;


D = 1/2*rho*v^2*A_nikha*Cd;

xprime = [Th_nikha/m-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-omega^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)); ...
-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+omega^2*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta)); ...
v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+omega^2*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta)); ...
v*sin(phi); ...
v/r*cos(phi)*cos(A); ...
v*cos(phi)*sin(A)/r*cos(delta); ...
mdot_nikha];
