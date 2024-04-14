
function Y = atmosphere(h,Lc)

global gamma R_air

R=R_air;		% ...Sea level gas constant for air [J/Kg-K]
g0=9.806;		% ...Sea level acceleration due to gravity [m/s^2]
Na=6.0220978e23;	% ...Avogadro's number
sigma=3.65e-10;		% ...Collision diameter [m] for air
S=110.4;		% ...Sutherland's Temperature
M0=28.964;		% ...Sea level molecular weight [g/mole]
T0=288.15;		% ...Sea level temperature [K]
P0=1.01325e5;		% ...Sea level pressure [N/m^2]
re=6378.14e3;		% ...Earth mean radius [m]
Beta=1.458e-6;		% ...Sutherland's constant [kg/m.s.K^0.5]

B=2/re;
layers=21;
Z=1e3*[0.00; 11.0191; 20.0631; 32.1619; 47.3501; 51.4125; 71.8020; 86; 100;
 110; 120; 150; 160; 170; 190; 230; 300; 400; 500; 600; 700; 2000];
T=[T0; 216.65; 216.65; 228.65; 270.65; 270.65; 214.65; 186.946;
 210.65; 260.65; 360.65; 960.65; 1110.60; 1210.65; 1350.65; 1550.65; 1830.65;
 2160.65; 2420.65; 2590.65; 2700.00; 2700.00];
M=[M0; 28.964; 28.964; 28.964; 28.964; 28.964; 28.962; 28.962; 28.880; 28.560;
 28.070; 26.920; 26.660; 26.500; 25.850; 24.690; 22.660; 19.940; 17.940; 16.840;
 16.170; 16.170];
a=[-6.5e-3; 0; 1e-3; 2.8e-3; 0; -2.8e-3; -2e-3;
 1.693e-3; 5.00e-3; 1e-2; 2e-2; 1.5e-2; 1e-2; 7e-3; 5e-3; 4e-3;
 3.3e-3; 2.6e-3; 1.7e-3; 1.1e-3; 0];

rho0=P0/(R*T0);
P(1)=P0;
T(1)=T0;
rho(1)=rho0;

for i=1:layers
  if ~(a(i) == 0)
	C1 = 1 + B*(T(i)/a(i) - Z(i));
	C2 = C1*g0/(R*a(i));
	C3 = T(i+1)/T(i);
	C4 = C3^(-C2);
	C5 = exp(g0*B*(Z(i+1)-Z(i))/(R*a(i)));
	C6 = C2+1;
	P(i+1)=P(i)*C4*C5;
	rho(i+1)=rho(i)*C5*C3^(-C6);
  else
	C7=-g0*(Z(i+1)-Z(i))*(1-B*(Z(i+1)+Z(i))/2)/(R*T(i));
	P(i+1)=P(i)*exp(C7);
	rho(i+1)=rho(i)*exp(C7);
  end
end

for i= 1:21
  if h< Z(i+1)
    if ~(a(i)==0)
	C1=1+B*(T(i)/a(i)-Z(i));
	TM=T(i)+a(i)*(h-Z(i));
	C2=C1*g0/(R*a(i));
	C3=TM/T(i);
	C4=C3^(-C2);
	C5=exp(B*g0*(h-Z(i))/(R*a(i)));
	C6=C2+1;
	PR=P(i)*C4*C5; %...Static Pressure [N/m^2]
	rhoE=C5*rho(i)*C3^(-C6); %...Density [kg/m^3]
    else
	TM=T(i);
	C7=-g0*(h-Z(i))*(1-(h+Z(1))*B/2)/(R*T(i));
	PR=P(i)*exp(C7);  %...Static Pressure [N/m^2]
	rhoE=rho(i)*exp(C7);  %...Density [kg/m^3]
    end
  MOL = M(i) + (M(i+1)-M(i))*(h - Z(i))/(Z(i+1) - Z(i));
  TM  = MOL*TM/M0; %...Kinetic Temperature
  asound = sqrt(gamma*R*TM); %...Speed of sound [m/s]
  Vm = sqrt(8*R*TM/pi);
  m = MOL*1e-3/Na;
  n = rhoE/m;
  F = sqrt(2)*pi*n*sigma^2*Vm;
  Lamda = Vm/F;  % Mean Free path [m]
  Kn = Lamda/Lc;  % Knudsen number
  Y = [rhoE; asound; Kn];
  return;
  end
end

