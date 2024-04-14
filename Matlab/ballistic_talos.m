function xprime = ballistic_talos(t,x);

global g0 rho0 hscale Re Ma15C subs trans a b c d e f A_talos


%This function computes the derivatives for a ballistic flight into the atmosphere.

Ma = x(1)/Ma15C;

if Ma < subs
	Cd = a;
    elseif Ma < trans
	Cd = -b+c*Ma;
    else
	Cd = d+e/sqrt(Ma^2-f);
end

g = g0/(1+x(4)/Re)^2;
rho = rho0*exp(-x(4)/hscale);

%Drag calculation

D = 1/2*rho*x(1)^2*A_talos*Cd;

xprime = [-D/x(7)-g*sin(x(3))-omega^2*x(4)*cos(x(5))*(cos(x(3))*cos(x(2))*sin(x(5))-sin(x(3))*cos(x(5))); ...
x(1)*cos(x(3))*sin(x(2))*tan(x(5))/x(4)+omega^2*x(4)*sin(x(2))*sin(x(5))*cos(x(5))/x(1)*cos(x(3))-2*omega/cos(x(3))*(sin(x(3))*cos(x(2))*cos(x(5))-cos(x(3))*sin(x(5))); ...
x(1)/x(4)*cos(x(3))-g*cos(x(3))/x(1)+2*omega*sin(x(2))*cos(x(5))+omega^2*x(4)*cos(x(5))/x(1)*(sin(x(3))*cos(x(2))*sin(x(5))+cos(x(3))*cos(x(5))); ...
x(1)*sin(x(3)); ...
x(1)/x(4)*cos(x(3))*cos(x(2)); ...
x(1)*cos(x(3))*sin(x(2))/x(4)*cos(x(5)); ...
0];
