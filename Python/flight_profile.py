#Debrisk UdeA module

from numpy import *
import numpy as np

from scipy.integrate import odeint

from math import *

from scipy import interpolate


diam_fustrum	= 0.8							# Fustrum diameter [m]
A_fustrum		= (pi/4)*(diam_fustrum)**2	# Fustrum area [m^2]
Lc_fustrum		= 0.2							
diam_capsule	= 0.78							# Reentry capsule diameter [m]
A_capsule	=(pi/4)*(diam_capsule)**2	


#%%
deg	= pi/180			# Convert degrees to radians
rho0	= 1.225				# Sea level atmosphere density
gamma	= 1.4				# Sea level specific heat ratio
R_air	= 257.058			# Sea level gas constant for air [J/Kg-K]

# NGA/NASA EGM96, N=M=360, WGS84 REVISED EARTH GRAVITATIONAL MODEL
Re	= 6378.137e3						# Semi-major axis of WGS 84 Ellipsoid (Radius of the earth equator for the NGA/NASA EGM96 WGS84) Gravitational model[m]
R_polar	= 6356.752314245e3		# Semi-minor axis of WGS 84 Ellipsoid (Polar radius of the earth) [m]
epsilon	= 1.0/298.2572235630	# Flattening of WGS 84 Ellipsoid
omega	= 7.2921150e-5				# Angular speed of Earth's rotation [rad/sec]
mu	= 3.986004418e14 				# mu=GMe
J2	= 1.08263e-3						# J2, J3, J4 harmonics for modelling the WGS84
J3	= 2.532153e-7
J4	= 1.6109876e-7
g = 9.80665 							# standard Earth gravity

# FLOW REGIME CONSTANTS FOR DRAG CALCULATION

subs	= 0.89				# Limit of subsonic regime for this rocket nose (Hoerner 1965)
trans	= 1.13				# High limit of transonic regime for this rocket nose (Hoerner 1965)

# Empirical coeficients for Drag calculation (from Hoerner, 1965)
a	= 0.1		
b	= 1
c	= 1.231
d	= 0.141
e	= 0.129
f	= 1

#%% 
def powered(x,t,q1,q2,q3): #  (A_talos,Th_talos,mdot_talos1) =(q1,q2,q3) tuple qi parameter 
                           #x = z1_0  (t,x,p)
    

    hturn = 100 # At this altitude begins the pitchover [m]

    v=x[0] # Parameters brought by vector x, are called with mnemotechnic names
    A=x[1]
    phi=x[2]
    r=x[3]
    delta=x[4]
    lamda=x[5]
    m=x[6]
    
    A_rocket=q1 # Rocket Parameters brought by vector p,are called with mnemotechnic names
    Th_rocket=q2
    mdot_rocket=q3
    
    R = Re*(1-epsilon*(sin(delta))**2)
    
    Lc = Lc_fustrum
    h = r-R
    Y = atmosphere(h,Lc)
    rho = Y[0]
    v_sound = Y[1]
    Knudsen = Y[2]
    Ma = v/v_sound
    s = Ma*sqrt(gamma/2)
    continuum_flow = 0.01
    molecular_flow = 10
    
    if Ma < subs:
        Cd_cont = a
    elif Ma < trans:
            Cd_cont = -b+c*Ma
    else:
        Cd_cont = d+e/sqrt(Ma**2-f)

    if Knudsen <= continuum_flow:
        Cd = Cd_cont    
    else:
        Cd_freemol = 1.75+sqrt(pi)/(2*s)
        if Knudsen > molecular_flow:
            Cd = Cd_freemol
        else:
            Cd = Cd_cont+(Cd_freemol-Cd_cont)*(0.333*log10(Knudsen/sin(pi/6))+0.5113)
        
    [g_center, g_delta]=gravity(r,delta)
    
    gc=abs(-g_center) # This action reverses the gravity vector sense.
    gdelta=-g_delta
    
    D = (1/2)*rho*(v**2)*A_fustrum*Cd
    
    if h<=hturn:
        xprime = [Th_rocket/m-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-(omega**2)*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)) \
                  ,0 
                  ,0  
                  ,v*sin(phi)
                  ,v/r*cos(phi)*cos(A)
                  ,v*cos(phi)*sin(A)/r*cos(delta)
                  ,mdot_rocket
                  ,-D/m 
                  ,-abs(-gc*sin(phi)/m)]
        
        
    
    else:
        xprime = [Th_rocket/m-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-(omega**2)*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)) \
                  ,-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+(omega**2)*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta))
                  ,v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+(omega**2)*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta))
                  ,v*sin(phi)
                  ,v/r*cos(phi)*cos(A)
                  ,v*cos(phi)*sin(A)/r*cos(delta)
                  ,mdot_rocket
                  ,-D/m
                  ,-abs(-gc*sin(phi)+gdelta*cos(phi)*cos(A))/m ]
        
    
        
        
    return xprime 
#%%
def atmosphere(h,Lc):
    

   
    R=R_air		# ...Sea level gas constant for air [J/Kg-K]
    g0=9.806		# ...Sea level acceleration due to gravity [m/s^2]
    Na=6.0220978e23	# ...Avogadro's number
    sigma=3.65e-10		# ...Collision diameter [m] for air
    S=110.4		# ...Sutherland's Temperature
    M0=28.964		# ...Sea level molecular weight [g/mole]
    T0=288.15		# ...Sea level temperature [K]
    P0=1.01325e5		# ...Sea level pressure [N/m^2]
    re=6378.14e3		# ...Earth mean radius [m]
    Beta=1.458e-6		# ...Sutherland's constant [kg/m.s.K^0.5]
    
    B=2/re
    layers=21
    
    Z= array([0.00, 11.0191, 20.0631, 32.1619, 47.3501, 51.4125, 71.8020, 86, 100,
     110, 120, 150, 160, 170, 190, 230, 300, 400, 500, 600, 700, 2000])*1e3
            
    T= array([T0, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946,
     210.65, 260.65, 360.65, 960.65, 1110.60, 1210.65, 1350.65, 1550.65, 1830.65,
     2160.65, 2420.65, 2590.65, 2700.00, 2700.00])
            
    M= array([M0, 28.964, 28.964, 28.964, 28.964, 28.964, 28.962, 28.962, 28.880, 28.560,
     28.070, 26.920, 26.660, 26.500, 25.850, 24.690, 22.660, 19.940, 17.940, 16.840,
     16.170, 16.170])
            
    a= array([-6.5e-3, 0, 1e-3, 2.8e-3, 0, -2.8e-3, -2e-3,
     1.693e-3, 5.00e-3, 1e-2, 2e-2, 1.5e-2, 1e-2, 7e-3, 5e-3, 4e-3,
     3.3e-3, 2.6e-3, 1.7e-3, 1.1e-3, 0])
    
    rho0=P0/(R*T0)
    P=[P0]
    #T(1)=T0;
    rho=[rho0]
    
    for i in range(layers):
        if (a[i] != 0):
            C1 = 1 + B*(T[i]/a[i] - Z[i])
            C2 = C1*g0/(R*a[i])
            C3 = T[i+1]/T[i]
            C4 = C3**(-C2)
            C5 = exp(g0*B*(Z[i+1]-Z[i])/(R*a[i]))
            C6 = C2+1
            P.append(P[len(P)-1]*C4*C5)
            rho.append(rho[len(rho)-1]*C5*C3**(-C6))
            #P(i+1)=P(i)*C4*C5
            #rho(i+1)=rho(i)*C5*C3^(-C6)
        else:
            C7 = -g0*(Z[i+1]-Z[i])*(1-B*(Z[i+1]+Z[i])/2)/(R*T[i])
            P.append(P[len(P)-1]*exp(C7))
            rho.append(rho[len(rho)-1]*exp(C7))
            
            #P(i+1)=P(i)*exp(C7)
            #rho(i+1)=rho(i)*exp(C7)
      
    for i in range(layers):
        if h< Z[i+1]:
            if (a[i]!=0):
                C1 = 1 + B*(T[i]/a[i] - Z[i])
                TM = T[i]+a[i]*(h-Z[i])
                C2 = C1*g0/(R*a[i])
                C3 = TM/T[i]
                C4 = C3**(-C2)
                C5 = exp(B*g0*(h-Z[i])/(R*a[i]))
                C6 = C2+1
                PR = P[i]*C4*C5 #...Static Pressure [N/m^2]
                rhoE = C5*rho[i]*C3**(-C6) #...Density [kg/m^3]
            else:
                TM=T[i]
                C7 = -g0*(h-Z[i])*(1-(h+Z[1])*B/2)/(R*T[i])
                PR = P[i]*exp(C7)  #...Static Pressure [N/m^2]
                rhoE = rho[i]*exp(C7)  #...Density [kg/m^3]
        
            MOL = M[i] + (M[i+1]-M[i])*(h - Z[i])/(Z[i+1] - Z[i])
            TM  = MOL*TM/M0 #...Kinetic Temperature
            asound = sqrt(gamma*R*TM) #...Speed of sound [m/s]
            Vm = sqrt(8*R*TM/pi)
            m = MOL*1e-3/Na
            n = rhoE/m
            F = sqrt(2)*pi*n*(sigma**2)*Vm
            Lamda = Vm/F  # Mean Free path [m]
            Kn = Lamda/Lc  # Knudsen number
            Y = [rhoE, asound, Kn]
            return Y
    

#%%
    
def gravity(r,latitude):
    
    phi=pi/2-latitude

    g_center=-mu*(1-1.5*J2*(3*(cos(phi)**2)-1)*(Re/r)**2-2*J3*cos(phi)*(5*cos(phi)**2-3)*
                  (Re/r)**3-(5/8)*J4*(35*cos(phi)**4 -30*cos(phi)**2+3)*(Re/r)**4)/r**2

    g_north=3*mu*sin(phi)*cos(phi)*(Re/r)*(Re/r) *(J2+0.5*J3*(5*cos(phi)**2-1)*
                     (Re/r)/cos(phi)+(5/6)*J4*(7*cos(phi)**2-1)*(Re/r)**2)/r**2
    
    return [g_center, g_north ]


#%%
    
def ballistic(x,t): #  (A_talos,Th_talos,mdot_talos1) =(q1,q2,q3) tuple qi parameter 
                           #x = z1_0  (t,x,p)

    v=x[0] # Parameters brought by vector x, are called with mnemotechnic names
    A=x[1]
    phi=x[2]
    r=x[3]
    delta=x[4]
    lamda=x[5]
    m=x[6]

    R = Re*(1-epsilon*(sin(delta))**2)
    
    Lc = Lc_fustrum
    h = r-R
    Y = atmosphere(h,Lc)
    rho = Y[0]
    v_sound = Y[1]
    Knudsen = Y[2]
    Ma = v/v_sound
    s = Ma*sqrt(gamma/2)
    continuum_flow = 0.01
    molecular_flow = 10
    
    if Ma < subs:
        Cd_cont = a
    elif Ma < trans:
            Cd_cont = -b+c*Ma
    else:
        Cd_cont = d+e/sqrt(Ma**2-f)
    
    if Knudsen <= continuum_flow:
        Cd = Cd_cont
    else:
        Cd_freemol = 1.75+sqrt(pi)/(2*s)
        if Knudsen > molecular_flow:
          # Modificacion 
        
            Cd = Cd_freemol
        else:
            Cd = Cd_cont+(Cd_freemol-Cd_cont)*(0.333*log10(Knudsen/sin(pi/6))+0.5113)
    
    
    [g_center, g_delta]=gravity(r,delta)
    
    gc=abs(-g_center) # This action reverses the gravity vector sense.
    gdelta=-g_delta
    
    D = (1/2)*rho*(v**2)*A_fustrum*Cd
    
    """
    xprime = [-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-(omega**2)*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta))
              ,-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+(omega**2)*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta))
              ,v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+(omega**2)*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta))
              ,0#v*sin(phi)
              ,0#v/r*cos(phi)*cos(A)
              ,0#v*cos(phi)*sin(A)/r*cos(delta)
              ,0
              ,0#-D/m
              ,0#-abs(-gc*sin(phi)+gdelta*cos(phi)*cos(A))/m]
              ]
    """
    xprime = [-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-(omega**2)*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)) \
              ,-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+(omega**2)*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta))
              ,v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+(omega**2)*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta))
              ,v*sin(phi)
              ,v/r*cos(phi)*cos(A)
              ,v*cos(phi)*sin(A)/r*cos(delta)
              ,0
              ,-D/m
              ,-abs(-gc*sin(phi)+gdelta*cos(phi)*cos(A))/m]
    
    
    return xprime 


def reentry(x,t):
    
    v=x[0] # Parameters brought by vector x, are called with mnemotechnic names
    A=x[1]
    phi=x[2]
    r=x[3]
    delta=x[4]
    lamda=x[5]
    m=x[6]

    R = Re*(1-epsilon*(sin(delta))**2)
    
    Lc = Lc_fustrum
    h = r-R
    Y = atmosphere(h,Lc)
    rho = Y[0]
    v_sound = Y[1]
    Knudsen = Y[2]
    Ma = v/v_sound
    s = Ma*sqrt(gamma/2)
    continuum_flow = 0.01
    molecular_flow = 10
    
    
    Cd_cont = Cd_reentry(Ma)
    Cd_freemol = 1.75+sqrt(pi)/(2*s)
        
    
    if Knudsen <= continuum_flow:
        Cd = Cd_cont
        
    elif Knudsen >= molecular_flow:
        Cd = Cd_freemol
    
    else:
        Cd = Cd_cont+(Cd_freemol-Cd_cont)*(0.333*log10(Knudsen/sin(pi/6))+0.5113)
    
       
    
    
    [g_center, g_delta]=gravity(r,delta)
    
    gc=abs(-g_center) # This action reverses the gravity vector sense.
    gdelta=-g_delta
    
    D = (1/2)*rho*(v**2)*A_capsule*Cd
    
    xprime = [-D/m-gc*sin(phi)+gdelta*cos(phi)*cos(A)-(omega**2)*r*cos(delta)*(cos(phi)*cos(A)*sin(delta)-sin(phi)*cos(delta)) \
              ,-gdelta*sin(A)/v*cos(phi)+v*cos(phi)*sin(A)*tan(delta)/r+(omega**2)*r*sin(A)*sin(delta)*cos(delta)/v*cos(phi)-2*omega*(tan(phi)*cos(A)*cos(delta)-sin(delta))
              ,v/r*cos(phi)-gc*cos(phi)/v-gdelta*sin(phi)*cos(A)/v+2*omega*sin(v)*cos(delta)+(omega**2)*r*cos(delta)/v*(sin(phi)*cos(A)*sin(delta)+cos(phi)*cos(delta))
              ,v*sin(phi)
              ,v/r*cos(phi)*cos(A)
              ,v*cos(phi)*sin(A)/r*cos(delta)
              ,0]
    
    return xprime 
    
    
def Cd_reentry(mach):
    

    mach_number = [0, 0.3, 0.5, 0.8, 0.81, 0.9, 0.96, 1.01, 1.2, 1.21, 1.6, 2, 3.5, 4.2, 6, 8, 99]
    Cd_number = [0.41, 0.41, 0.41 ,0.45, 0.49 ,0.5 ,0.58, 0.62, 0.71, 0.71 ,0.68 ,0.66, 0.6587 ,0.6482, 0.6275, 0.6, 0.6]
    Cd_aux = interpolate.interp1d(mach_number,Cd_number)
    Cd = Cd_aux(mach)
    return Cd
    

#%%
 #function s = distVincenty(lat1,lat2,lon1,lon2)
 
def distVincenty(lat1,lat2,lon1,lon2):
    

    a=Re
    b=R_polar  # NGA/NASA EGM96, N=M=360, WGS84 Parameters
    L=lon2-lon1
    U1=atan((1-epsilon)*tan(lat1))  # Reduced latitude
    U2=atan((1-epsilon)*tan(lat2))
    sinU1=sin(U1)
    sinU2=sin(U2)
    cosU1=cos(U1)
    cosU2=cos(U2)
    
    lamda=L 
    lamdaP=0 
    iterlimit=100
    
    while (abs(lamda-lamdaP)>1e-12 and iterlimit>0):
        sinlamda=sin(lamda)
        coslamda=cos(lamda)
        sinsigma=sqrt((cosU2*sinlamda)**2+(cosU1*sinU2-sinU1*cosU2*coslamda)**2)
        if (sinsigma==0): 
            s=0 
            return s
    	# Coincident points
        cossigma=sinU1*sinU2+cosU1*cosU2*coslamda
        sigma=atan2(sinsigma,cossigma)
        sinalpha=cosU1*cosU2*sinlamda/sinsigma
        cossqalpha=1-sinalpha**2
        cos2sigmaM=cossigma-2*sinU1*sinU2/cossqalpha
        if (isnan(cos2sigmaM)==True):
            cos2sigmaM=0 
    	
        C=epsilon/16*cossqalpha*(4+epsilon*(4-3*cossqalpha))
        lamdaP=lamda
        lamda=L+(1-C)*epsilon*sinalpha*(sigma+C*sinsigma*(cos2sigmaM+C*cossigma*(-1+2*cos2sigmaM**2)))
        if (iterlimit==0):
            s=nan
            return s

        iterlimit -=1
        
    
    
    uSq=cossqalpha*(a**2-b**2)/b**2
    A=1+uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B=uSq/1024*(256+uSq*(-128+uSq*(74-47*uSq)))
    deltasigma=B*sinsigma*(cos2sigmaM+B/4*(cossigma*(-1+2*cos2sigmaM**2)-(B/6)*cos2sigmaM*(-3+4*sinsigma**2)*(-3+4*cos2sigmaM**2)))
    s=b*A*(sigma-deltasigma)
    return s 






#%%
def ned2ecefv(uNorth, vEast,wDown, lat0, lon0):
    wUp = wDown *-1
    cosPhi = np.cos(lat0)
    sinPhi = np. sin(lat0)
    cosLambda = np.cos(lon0)
    sinLambda = np.sin(lon0)

    t = cosPhi*wUp - sinPhi*uNorth
    w = sinPhi*wUp + cosPhi*uNorth

    u = cosLambda*t - sinLambda*vEast
    v = sinLambda*t + cosLambda*vEast
    return [u,v,w] # u = Vx , v = Vy, w = Vz

#%%

def dynamic_pressure(v,h):
    

# This function takes the velocity and altitude... and computes
# the dynamic pressure during ascent and reentry [kPa]

    H= array(h)*1.e3  # reverses the altitude from km to m.
    V= array(v)*1.e3
    Lc=1 # Input argument for atmosphere function, it doesn't matter here.
    q =[]
    for i in range(len(V)):
        Y = atmosphere(H[i],Lc)
        #Modificar 
        Rho=Y[0]
        q.append(1/2*Rho*V[i]**2)
      	
    
    Q=(array(q)/1000)
    return Q


       
def ballistic_event(z_0,t0):
    
    i=0
    for dt in arange(100,10000,100):
        tf=t0+dt
        tspan =linspace(t0,tf,300)       
        z_t=odeint(ballistic,z_0,tspan) # Options 1 
        if i==0:
            z=copy(z_t)
            
        else:
            cnt = 0
            for i in (z_t[:,2]):
                if i <0:
                    z=concatenate((z,z_t[0:cnt+1]))
                    return [z,tspan[len(z_t[0:cnt+1])]]
                cnt +=1
            z=concatenate((z,z_t))
        t0=tf
        z_0=z_t[-1,:]
        i+=1
 
def reentry_event(z_0,t0):
    
    i=0
    for dt in arange(100,10000,100):
        tf=t0+dt
        tspan =linspace(t0,tf,300)       
        z_t=odeint(reentry,z_0,tspan)
        R=Re*(1-epsilon*(sin(z_t[-1,4]))**2)+5000 
        if i==0:
            z=copy(z_t)
            
        else:
            cnt = 0
            for i in (z_t[:,3]):
                if i <R:
                    z=concatenate((z,z_t[0:cnt+1]))
                    return [z,tspan[len(z_t[0:cnt+1])]]
                cnt +=1
            z=concatenate((z,z_t))
        t0=tf
        z_0=z_t[-1,:]
        i+=1         


