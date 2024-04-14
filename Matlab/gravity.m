function [g_center g_north]=gravity(r,latitude)

global Re mu J2 J3 J4
% Geophysical constants for the NGA/NASA EGM96 Gravitational model
phi=pi/2-latitude;

g_center=-mu*(1-1.5*J2*(3*cos(phi)^2-1)*(Re/r)^2-2*J3*cos(phi)*(5*cos(phi)^2-3)*(Re/r)^3-(5/8)*J4*(35*cos(phi)^4 ...
-30*cos(phi)^2+3)*(Re/r)^4)/r^2;

g_north=3*mu*sin(phi)*cos(phi)*(Re/r)*(Re/r) ...
*(J2+0.5*J3*(5*cos(phi)^2-1) ...
*(Re/r)/cos(phi) ...
+(5/6)*J4*(7*cos(phi)^2-1)*(Re/r)^2)/r^2;

