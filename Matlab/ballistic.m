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

if Ma < subs
   Cd_cont = a;
   elseif Ma < trans
   Cd_cont = -b+c*Ma;
   else
   Cd_cont = d+e/sqrt(Ma^2-f);
end
Cd_freemol = 1.75+sqrt(pi)/(2*s);

if Knudsen <= continuum_flow
   Cd = Cd_cont; 
elseif Knudsen >= molecular_flow
   Cd = Cd_freemol;
else
   Cd = Cd_cont+(Cd_freemol-Cd_cont)*(0.333*log10(Knudsen/sin(pi/6))+0.5113);
end

[g_center g_delta]=gravity(r,delta);

gc=abs(-g_center); % This action reverses the gravity vector sense.
gdelta=-g_delta;


D = 1/2*rho*v^2*A_fustrum*Cd;

xprime = [-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-omega^2*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)); ...
-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+omega^2*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta)); ...
v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+omega^2*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta)); ...
v*sin(phi); ...
v/r*cos(phi)*cos(A); ...
v*cos(phi)*sin(A)/r*cos(delta); ...
0; ...
-D/m; ...
-abs(-gc*sin(phi)+gdelta*cos(phi)*cos(A))/m];
