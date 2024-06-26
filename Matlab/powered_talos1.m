function xprime = powered_talos1(t,x);

global g0 rho0 hscale Re Ma15C subs trans a b c d e f A_talos Th_talos mdot_talos omega epsilon

h = x(4)-Re;

Y = atmosphere(h);

rho = Y(1);

Ma = x(1)/Y(2);

if Ma < subs
   Cd = a;
   elseif Ma < trans
	Cd = -b+c*Ma;
    else
	Cd = d+e/sqrt(Ma^2-f);
end

[g_center g_delta]=gravity(x(4),x(5));

gc=abs(-g_center); % This action reverses the gravity vector sense.
gdelta=-g_delta;

D = 1/2*rho*x(1)^2*A_talos*Cd;

xprime = [Th_talos/x(7)-D/x(7)-gc*sin(x(3))+gdelta*cos(x(3))*cos(x(2))-omega^2*x(4)*cos(x(5))*(cos(x(3))*cos(x(2))*sin(x(5))-sin(x(3))*cos(x(5))); ...
0; ...
0; ...
x(1)*sin(x(3)); ...
x(1)/x(4)*cos(x(3))*cos(x(2)); ...
x(1)*cos(x(3))*sin(x(2))/x(4)*cos(x(5)); ...
mdot_talos];
