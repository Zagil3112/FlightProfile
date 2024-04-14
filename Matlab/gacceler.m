function Gs = gacceler(p,v,A,phi,r,delta,m)

global gamma Re subs trans a b c d e f omega epsilon

Th_rocket=p(1);
Area=p(2);
Lc=p(3);
flightphase=p(4);

for i=1:length(A)
	R=Re*(1-epsilon*(sin(delta(i)))^2);
	h=r(i)-R;
	Y=atmosphere(h,Lc);
	rho=Y(1);
	vsound=Y(2);
	knudsen=Y(3);
	Ma=v(i)/vsound;
	s=Ma*sqrt(gamma/2);
	continuum_flow=0.01;
	molecular_flow=10;
	if flightphase==1  % ascent phase
		if Ma < subs
   		Cd_cont = a;
   	elseif Ma < trans
   	Cd_cont = -b+c*Ma;
   	else
   	Cd_cont = d+e/sqrt(Ma^2-f); 
		end
	else
	Cd_cont = Cd_reentry1221(Ma); % reentry phase
	end
	Cd_freemol = 1.75+sqrt(pi)/(2*s);
	if knudsen <= continuum_flow
   	Cd = Cd_cont; 
		elseif knudsen >= molecular_flow
   	Cd = Cd_freemol;
		else
   	Cd = Cd_cont+(Cd_freemol-Cd_cont)*(0.333*log10(knudsen/sin(pi/6))+0.5113);
	end
	
	[g_center g_delta]=gravity(r(i),delta(i));
	gc=abs(-g_center); % This action reverses the gravity vector sense.
	gdelta=-g_delta;
	
	D = 1/2*rho*v(i)^2*Area*Cd;

	Gxv(i)=Th_rocket/m(i)-D/m(i)-gc*sin(phi(i))+gdelta*cos(phi(i))*cos(A(i));
	Gyv(i)=gdelta*sin(A(i));
	Gzv(i)=gc*cos(phi(i))+gdelta*sin(phi(i))*cos(A(i));
	Gtotalsq=Gxv(i)^2+Gyv(i)^2+Gzv(i)^2;
	Gtotal(i)=sqrt(Gtotalsq);
end

Gs=[Gxv' Gyv' Gzv' Gtotal'];
