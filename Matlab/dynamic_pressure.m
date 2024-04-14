function Q=dynamic_pressure(v,h)

% This function takes the velocity and altitude... and computes
% the dynamic pressure during ascent and reentry [kPa]

H=h*1.e3;  % reverses the altitude from km to m.
V=v*1.e3;
Lc=1; % Input argument for atmosphere function, it doesn't matter here.

  for i=1:length(v)
      	Y=atmosphere(H(i),Lc);
	Rho=Y(1);
	q(i)=1/2*Rho*V(i)^2;
  end	

Q=(q.')/1000;
