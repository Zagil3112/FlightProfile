function s = distVincenty(lat1,lat2,lon1,lon2)

global Re R_polar epsilon

a=Re; b=R_polar;  % NGA/NASA EGM96, N=M=360, WGS84 Parameters
L=lon2-lon1;
U1=atan((1-epsilon)*tan(lat1));  % Reduced latitude
U2=atan((1-epsilon)*tan(lat2));
sinU1=sin(U1);
sinU2=sin(U2);
cosU1=cos(U1);
cosU2=cos(U2);

lamda=L; lamdaP=0; iterlimit=100;

while (abs(lamda-lamdaP)>1e-12 & --iterlimit>0)
	sinlamda=sin(lamda); coslamda=cos(lamda);
	sinsigma=sqrt((cosU2*sinlamda)^2+(cosU1*sinU2-sinU1*cosU2*coslamda)^2);
	if (sinsigma==0) 
	  s=0; 
	  return;
	end % Coincident points
	cossigma=sinU1*sinU2+cosU1*cosU2*coslamda;
	sigma=atan2(sinsigma,cossigma);
	sinalpha=cosU1*cosU2*sinlamda/sinsigma;
	cossqalpha=1-sinalpha^2;
	cos2sigmaM=cossigma-2*sinU1*sinU2/cossqalpha;
	if (isnan(cos2sigmaM)==1) 
	  cos2sigmaM=0; 
	end
	C=epsilon/16*cossqalpha*(4+epsilon*(4-3*cossqalpha));
	lamdaP=lamda;
	lamda=L+(1-C)*epsilon*sinalpha*(sigma+C*sinsigma*(cos2sigmaM+C*cossigma*(-1+2*cos2sigmaM^2)));
end

if (iterlimit==0)
  s=NaN;
  return;
end


uSq=cossqalpha*(a^2-b^2)/b^2;
A=1+uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
B=uSq/1024*(256+uSq*(-128+uSq*(74-47*uSq)));
deltasigma=B*sinsigma*(cos2sigmaM+B/4*(cossigma*(-1+2*cos2sigmaM^2)-B/6 ...
*cos2sigmaM*(-3+4*sinsigma^2)*(-3+4*cos2sigmaM^2)));
s=b*A*(sigma-deltasigma);


