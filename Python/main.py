# -*- coding: utf-8 -*-
import numpy as np
from math import *  
from flight_profile import *
from pandas import DataFrame
from scipy.integrate import odeint
import matplotlib.pyplot as mt 
# ...THIS PROGRAM CALCULATES THE TRAJECTORY FOR A THREE STAGE SOUNDING ROCKET

from sys import exit


# DECLARE GLOBAL VARIABLES FOR USE IN EVERY FUNCTION
#global rho0 gamma R_air Re R_polar subs trans a b c d e f omega diam_talos Lc_capsule A_capsule Lc_fustrum A_fustrum epsilon mu J2 J3 J4

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

# ROCKET CONSTANTS FOR TALOS-1 STAGE
diam_talos	= 0.76									# Vehicle Diameter [m]
A_talos		= (pi/4)*(diam_talos)**2				# Frontal area [m^2]
m0_talos	= 1996										# Talos total mass [kg]
mfinal_talos	= 496								# Talos Unfueled mass [kg]
Th_talos 	= 516e3									# Talos thrust [N]
mdot_talos	=-288.46								# Talos propellant consumption rate [kg/s]
mprop_talos	= m0_talos - mfinal_talos			# Talos Propellant mass [kg]
tburn_talos	= mprop_talos/abs(mdot_talos)	# Talos burning time [s]

# ROCKET CONSTANTS FOR TOMAHAWK STAGE
diam_toma	= 0.58									# Vehicle Diameter [m]
A_toma		= (pi/4)*(diam_toma)**2				# Frontal area [m^2]
m0_toma		= 1363									# Tomahawk total mass [kg]
mfinal_toma	= 602									# Tomahawk unfueled mass
Th_toma		= 457e3									# Tomahawk thrust [N]
mdot_toma	= -217.43								# Tomahawk propellant consumption rate [kg/s]
mprop_toma	= m0_toma - mfinal_toma			# Tomahawk propellant mass [kg]
tburn_toma	= mprop_toma/abs(mdot_toma)		# Tomahawk burning time [s]

# ROCKET CONSTANTS FOR NIKHA STAGE
diam_nikha		= 0.44									# Vehicle Diameter [m]
A_nikha			= (pi/4)*(diam_nikha)**2				# Frontal area [m^2]
m0_nikha			= 399									# Nikha total mass [kg]
mfinal_nikha	= 70										# Nikha unfueled mass [kg]
Th_nikha			= 50500									# Nikha thrust [N]
mdot_nikha		= -19.35								# Nikha propellant consumption rate [Kg/s]
mprop_nikha		= m0_nikha - mfinal_nikha			# Nikha propellant mass [kg]
tburn_nikha		= mprop_nikha/abs(mdot_nikha)	# Nikha burning time [s]

# CONSTANTS FOR ROCKET FUSTRUM
diam_fustrum	= 0.8							# Fustrum diameter [m]
A_fustrum		= (pi/4)*(diam_fustrum)**2	# Fustrum area [m^2]
Lc_fustrum		= 0.2							# Characteristic longitude for ascent Knudsen [m]

# CONSTANTS FOR CAPSULE
diam_capsule	= 0.78							# Reentry capsule diameter [m]
A_capsule	=(pi/4)*(diam_capsule)**2		# Reentry capsule frontal area [m]
Lc_capsule	= 0.21								# Characteristic longitude for reentry Knudsen [m]

# PAYLOAD MASS
m_pay		= 410 # [kg]

# INITIAL CONDITIONS
t0 = 0					# The flight time clock begins
tf1 = tburn_talos	# Talos burning time for time interval
#tspan1 = [t0,tf1]	# Time interval for first flight phase integration 
tspan1 = np.linspace(t0,tf1, 73)
v1st_0		= 0												# Initial velocity
A1st_0		= 85*deg										# Initial azimuth
phi1st_0	= 77*deg											# Initial angle wrt local horizontal
# prompt = 'Launch site latitude(deg): '
# lat = input(prompt)
# prompt2 = 'Launch site longitude(deg): '
# l_ong = input(prompt2)
# delta1st_0 = lat*deg
# lamda1st_0 = l_ong*deg

delta1st_0	= 37.940194*deg								# Launch site latitude
lamda1st_0	= -75.466389*deg								# Launch site longitude
r1st_0		= Re*(1-epsilon*(np.sin(delta1st_0))**2)	# Initial datum sea level due to WGS84 model
m1st_0		= m0_talos+m0_toma+m0_nikha+m_pay		# Total initial mass
vlossdrag1st_0	= 0											# Initial velocity loss	due to drag
vlossgrav1st_0	= 0											# Initial velocity loss due to gravity

z1_0=[v1st_0 ,A1st_0 ,phi1st_0, r1st_0 ,delta1st_0 ,lamda1st_0 ,m1st_0 ,vlossdrag1st_0 ,vlossgrav1st_0]

#%% Tolerancias 


# Specify tolerance and event stop at apogee for ode45 integrator
"""options1 = odeset('Events',@ballis_event,'reltol',1e-12,'abstol',1e-12)

# Specify tolerance and event stop at the end of the flight

options2 = odeset('Events',@reentry_event,'reltol',1e-12,'abstol',1e-12)"""

#%%
p1=(A_talos,Th_talos,mdot_talos) # Talos rocket parameters to be sent to debrisk.powered

z1 = odeint(powered,z1_0,tspan1, args =(p1))

#[t1,z1]=ode45(@powered,tspan1,z1_0,[],p1)

# Solution z(t) returned on the time interval [t0 tf] for first flight phase:
v1st				= z1[:,0]	 # Velocity [m/s]
A1st				= z1[:,1]	 # Azimuth angle [rad]
phi1st			   = z1[:,2]	 # Flight path angle [rad]
r1st				= z1[:,3]	 # Distance from Earth center [m]
delta1st			= z1[:,4]	 # Latitude [rad]
lamda1st			= z1[:,5]	 # Longitude [rad]
m1st				= z1[:,6]	 # Mass [kg]
vlossdrag1st	= z1[:,7]	 # Velocity loss due to drag [m/s]
vlossgrav1st	= z1[:,8]	 # Velocity loss due to gravity [m/s]

# Initial conditions for the second flight phase (1st ballistic)
# the same quantities as in  flight phase 1 apply

v2nd_0			= v1st[-1]			
A2nd_0			= A1st[-1]				
phi2nd_0			= phi1st[-1]					
r2nd_0			= r1st[-1]		
delta2nd_0		= delta1st[-1]					
lamda2nd_0		= lamda1st[-1]
m2nd_0			= m0_toma+m0_nikha+m_pay 	# Talos stage has been detached
vlossdrag2nd_0	= vlossdrag1st[-1]
vlossgrav2nd_0	= vlossgrav1st[-1]

t02	 = tf1			# Initial time for 2nd phase = final time 1st phase
tf2	 = t02+2		# Two seconds between the stages	
tspan2x = [t02, tf2]	# Time interval for second flight phase integration
tspan2 = np.linspace(t02,tf2,42)

z2_0 		= [v2nd_0, A2nd_0, phi2nd_0, r2nd_0, delta2nd_0 ,lamda2nd_0, m2nd_0 ,vlossdrag2nd_0, vlossgrav2nd_0]

z2 =odeint(ballistic,z2_0,tspan2)

#[t2,z2]=ode45(@ballistic,tspan2,z2_0)

# Solution z2(t) returned on the time interval [t02 tf2]
v2nd	   	 = z2[:,0]	# Velocity [m/s]
A2nd  		 = z2[:,1]	# Azimuth angle [rad]
phi2nd	    = z2[:,2]	# Flight path angle [rad]
r2nd	   	 = z2[:,3]# Distance from Earth Center [m]
delta2nd		 = z2[:,4]	# Latitude [rad]
lamda2nd		 = z2[:,5]	# Longitude [rad]
m2nd			 = z2[:,6]	# Mass [kg]
vlossdrag2nd = z2[:,7]	# Velocity loss due to drag [m/s]
vlossgrav2nd = z2[:,8]	# Velocity loss due to gravity [m/s]

# Initial conditions for the third flight phase (powered)
# # the same quantities as in  flight phase 1 and 2 apply

v3rd_0			= v2nd[-1]				
A3rd_0			= A2nd[-1]				
phi3rd_0			= phi2nd[-1]					
r3rd_0			= r2nd[-1]					
delta3rd_0		= delta2nd[-1]					
lamda3rd_0		= lamda2nd[-1]		
m3rd_0			= m2nd[-1]		
vlossdrag3rd_0	= vlossdrag2nd[-1]		
vlossgrav3rd_0	= vlossgrav2nd[-1]		

t03		= tf2  # Initial time for 3rd phase = 2nd phase final time
tf3		= t03+tburn_toma
tspan3x	= [t03, tf3]
tspan3 = np.linspace(t03,tf3,42)

z3_0		= [v3rd_0,A3rd_0 ,phi3rd_0 ,r3rd_0 ,delta3rd_0 ,lamda3rd_0 ,m3rd_0 ,vlossdrag3rd_0 ,vlossgrav3rd_0]

p2=(A_toma, Th_toma, mdot_toma) # Tomahawk rocket parameters to be sent to powered.m

z3 = odeint(powered,z3_0,tspan3, args =(p2))
#[t3,z3]=ode45(@powered,tspan3,z3_0,[],p2)

# Solution z3(t) returned on the time interval [t03 tf3]
# The same quantities as for previous phases apply

v3rd	   		= z3[:,0]
A3rd  			= z3[:,1]	
phi3rd	   	= z3[:,2]						
r3rd	   		= z3[:,3]						
delta3rd			= z3[:,4]						
lamda3rd			= z3[:,5]
m3rd				= z3[:,6]
vlossdrag3rd	= z3[:,7]
vlossgrav3rd	= z3[:,8]

# Initial conditions for the fourth flight phase (ballistic)
# The same quantities as for previous phases apply

v4th_0			= v3rd[-1]				
A4th_0			= A3rd[-1]				
phi4th_0			= phi3rd[-1]					
r4th_0			= r3rd[-1]					
delta4th_0		= delta3rd[-1]					
lamda4th_0		= lamda3rd[-1]
m4th_0			= m0_nikha+m_pay
vlossdrag4th_0	= vlossdrag3rd[-1]
vlossgrav4th_0	= vlossgrav3rd[-1]

t04	= tf3	# Initial time for 4th phase = 3rd phase final time
tf4	= t04+2 # Two seconds between the stages
tspan4X	= [t04,tf4]
tspan4 = np.linspace(t04,tf4,42)

z4_0		= [v4th_0, A4th_0, phi4th_0, r4th_0, delta4th_0 ,lamda4th_0 ,m4th_0 ,vlossdrag4th_0 ,vlossgrav4th_0]

z4=odeint(ballistic,z4_0,tspan4)
#[t4,z4]=ode45(@ballistic,tspan4,z4_0)

v4th	   	 = z4[:,0]						
A4th  		 = z4[:,1]						
phi4th	    = z4[:,2]						
r4th	   	 = z4[:,3]						
delta4th		 = z4[:,4]						
lamda4th		 = z4[:,5]
m4th			 = z4[:,6]
vlossdrag4th = z4[:,7]
vlossgrav4th = z4[:,8]

# Initial conditions for the fifth flight phase (Nikha-powered)
# The same quantities as for previous phases apply

v5th_0		   = v4th[-1]					
A5th_0		   = A4th[-1]				
phi5th_0		   = phi4th[-1]					
r5th_0		   = r4th[-1]					
delta5th_0	   = delta4th[-1]				
lamda5th_0	   = lamda4th[-1]
m5th_0		   = m4th[-1]
vlossdrag5th_0	= vlossdrag4th[-1]
vlossgrav5th_0	= vlossgrav4th[-1]

t05	 = tf4
tf5	 = t05+tburn_nikha
tspan5x= [t05 ,tf5]
tspan5 = np.linspace(t05,tf5,42)


z5_0	 = [v5th_0,A5th_0, phi5th_0 ,r5th_0 ,delta5th_0, lamda5th_0 ,m5th_0 ,vlossdrag5th_0 ,vlossgrav5th_0]

p3=(A_nikha, Th_nikha, mdot_nikha)
 # Nikha rocket parameters to be sent to powered.m function
z5=odeint(powered,z5_0,tspan5,args =(p3))
#[t5,z5]=ode45(@powered,tspan5,z5_0,[],p3)

v5th	   	 = z5[:,0]						
A5th  		 = z5[:,1]						
phi5th	    = z5[:,2]						
r5th	   	 = z5[:,3]						
delta5th		 = z5[:,4]						
lamda5th		 = z5[:,5]
m5th			 = z5[:,6]
vlossdrag5th = z5[:,7]
vlossgrav5th = z5[:,8]



# Initial conditions for the sixth flight phase (3rd ballistic)
# The same quantities as for previous phases apply

v6th_0			= v5th[-1]					
A6th_0			= A5th[-1]				
phi6th_0			= phi5th[-1]					
r6th_0			= r5th[-1]					
delta6th_0		= delta5th[-1]					
lamda6th_0		= lamda5th[-1]
m6th_0			= m_pay  # The nikha stage and the fustrum has been detached, only the capsule keeps on flight
vlossdrag6th_0	= vlossdrag5th[-1]
vlossgrav6th_0	= vlossgrav5th[-1]

t06	= tf5
tf6	= t06+10000    # Enough time to allow the capsule reachs the apogee point. 
							 # The integration will be stopped when ode45 reachs the event
							 # stop at apogee written at ballis_event.m
tspan6x= [t06, tf6]
#tspan6 = np.linspace(t06,tf6,100)
#tspan6 = np.linspace(t06,tf6,20973)

z6_0		= [v6th_0, A6th_0, phi6th_0, r6th_0, delta6th_0, lamda6th_0, m6th_0, vlossdrag6th_0, vlossgrav6th_0]


z6_aux = ballistic_event(z6_0,t06) # [z6,tf6]
z6 = z6_aux[0]
tspan6 = np.linspace(t06,z6_aux[1],len(z6))

#i=0
#for dt in np.arange(100,10000,50):
#    tf6=t06+dt
#    tspan6 =np.linspace(t06,tf6,100)
#    print("Integrando T0 =%F -- TF = %F"%(t06,tf6))
#    print("Condiciones iniciales")
#    print("")
#    print(z6_0)
#    print("")
#    z6_t=odeint(ballistic,z6_0,tspan6) # Options 1 
#    print("Condiciones finales ")
#    print(z6_t[-1,:])
#    
#    
#    t06=tf6
#    z6_0=z6_t[-1,:]
#    for i in (z6_t[:,2]):
#        if i <0:
#            
#        
#    if i==0:
#        z6=np.copy(z6_t)
#        tspan6_def =np.copy(tspan6)
#    else:
#        z6=np.concatenate((z6,z6_t))
#        tspan6_def =np.concatenate((tspan6_def,tspan6))
#    i+=1


#print(z6[-1])









































#i=0
#for dt in np.arange(100,10000,50):
#    tf6=t06+dt
#    tspan6 =np.linspace(t06,tf6,100)
#    print("Integrando T0 =%F -- TF = %F"%(t06,tf6))
#    print("Condiciones iniciales")
#    print("")
#    print(z6_0)
#    print("")
#    z6_t=odeint(ballistic,z6_0,tspan6) # Options 1 
#    print("Condiciones finales ")
#    print(z6_t[-1,:])
#    
#    
#    t06=tf6
#    z6_0=z6_t[-1,:]
#    if z6_t[-1,2]<-0.9:
#        break
#    if i==0:
#        z6=np.copy(z6_t)
#        tspan6_def =np.copy(tspan6)
#    else:
#        z6=np.concatenate((z6,z6_t))
#        tpan6_def =np.concatenate((tspan6_def,tspan6))
#    i+=1
#    

   


#print(tspan6)
#exit(0)

#i=0
#for dt in np.arange(100,10000,100):
#    tf6=t06+dt
#    tspan6 =np.linspace(t06,tf6,100)
#    print("Integrando ",t06,tf6)
#    print(z6_0)
#    z6_t=odeint(ballistic,z6_0,tspan6) # Options 1 
#    print(z6_t[-1,:])
#    R=Re*(1-epsilon*(np.sin(z6_t[-1,4]))**2)+5000;
#    print(R,z6_t[-1,3])
#    if z6_t[-1,3]<R:break
#    t06=tf6
#    z6_0=z6_t[-1,:]
#    if i==0:
#        z6=np.copy(z6_t)
#    else:
#        z6=np.concatenate((z6,z6_t))
#    i+=1
#    if i>1:break
#

#z6=odeint(ballistic,z6_0,tspan6)
#[t6,z6]=ode45(@ballistic,tspan6,z6_0,options1)

v6th	       = z6[:,0]						
A6th  	    = z6[:,1]						
phi6th	    = z6[:,2]						
r6th	       = z6[:,3]						
delta6th	    = z6[:,4]						
lamda6th	    = z6[:,5]
m6th		    = z6[:,6]
vlossdrag6th = z6[:,7]
vlossgrav6th = z6[:,8]

# Initial conditions for the seventh flight phase (reentry)
# The same quantities as for previous phases apply

v7th_0		= v6th[-1]
A7th_0		= A6th[-1]
phi7th_0		= phi6th[-1]
r7th_0		= r6th[-1]
delta7th_0	= delta6th[-1]
lamda7th_0	= lamda6th[-1]
m7th_0		= m6th[-1]

t07	 = tspan6[-1]
tf7	 = t07+10000   # Enough time to allow the capsule reachs the final flight point. 
							 # The integration will be stopped when ode45 reachs the event
							 # stop at 5 km. altitude over sea level, written at reentry_event.m
tspan7x = [t07, tf7]
#tspan7 = np.linspace(t07,tf7,34705)


z7_0 = [v7th_0 ,A7th_0 ,phi7th_0, r7th_0, delta7th_0, lamda7th_0 ,m7th_0]

z7_aux = reentry_event(z7_0,t07)

z7=z7_aux[0]
tspan7 = np.linspace(t07,z7_aux[1],len(z7))



#[t7,z7]=ode45(@reentry,tspan7,z7_0,options2)
#z7=np.zeros((10,7))

v7th	   = z7[:,0]						
A7th  	= z7[:,1]						
phi7th	= z7[:,2]						
r7th	   = z7[:,3]						
delta7th	= z7[:,4]						
lamda7th	= z7[:,5]
m7th		= z7[:,6]

# TOTAL FLIGHT DATA VECTORS

t 	 		 = 	np.concatenate((tspan1 ,tspan2, tspan3, tspan4, tspan5, tspan6, tspan7), axis =None)
v 	 		 = np.concatenate(( v1st ,v2nd, v3rd, v4th, v5th, v6th, v7th),axis =None)*1.e-3
v_asc		 = np.concatenate(( v1st, v2nd, v3rd, v4th ,v5th ,v6th),axis =None)*1.e-3
v_des		 = np.concatenate((	v7th),axis =None)*1.e-3
A	 		 = np.concatenate(( 	A1st ,A2nd, A3rd, A4th, A5th, A6th, A7th),axis =None)/deg
phi 	 	 = np.concatenate((	phi1st ,phi2nd, phi3rd, phi4th, phi5th, phi6th ,phi7th),axis =None)/deg
r	 		 = np.concatenate((	r1st ,r2nd ,r3rd ,r4th ,r5th, r6th, r7th),axis =None)*1.e-3
delta 	 = np.concatenate(( 	delta1st ,delta2nd, delta3rd, delta4th, delta5th, delta6th, delta7th),axis =None)/deg
del_aux = np.concatenate((		delta1st, delta2nd ,delta3rd ,delta4th ,delta5th, delta6th, delta7th),axis =None)
lamda	 	 = np.concatenate((	lamda1st, lamda2nd ,lamda3rd, lamda4th ,lamda5th, lamda6th ,lamda7th),axis =None)/deg
m 	 		 = np.concatenate(( 	m1st ,m2nd, m3rd, m4th, m5th, m6th, m7th),axis =None)
vlossdrag = np.concatenate((	vlossdrag1st, vlossdrag2nd ,vlossdrag3rd, vlossdrag4th ,vlossdrag5th ,vlossdrag6th),axis =None)*1.e-3
vlossgrav = np.concatenate((	vlossgrav1st ,vlossgrav2nd ,vlossgrav3rd, vlossgrav4th ,vlossgrav5th, vlossgrav6th),axis =None)*1.e-3

R1st	 	 =		Re*(1-epsilon*(np.sin((delta1st))**2))
R2nd	 	 =		Re*(1-epsilon*(np.sin(delta2nd))**2)
R3rd	 	 =		Re*(1-epsilon*(np.sin(delta3rd))**2)
R4th	 	 =		Re*(1-epsilon*(np.sin(delta4th))**2)
R5th	 	 =		Re*(1-epsilon*(np.sin(delta5th))**2)
R6th	 	 =		Re*(1-epsilon*(np.sin(delta6th))**2)
R7th	 	 =		Re*(1-epsilon*(np.sin(delta7th))**2)

h1st	 	 =		abs((r1st-R1st)*1.e-3)
h2nd	 	 =		abs((r2nd-R2nd)*1.e-3)
h3rd	 	 =		abs((r3rd-R3rd)*1.e-3)
h4th	 	 =		abs((r4th-R4th)*1.e-3)
h5th	 	 =		abs((r5th-R5th)*1.e-3)
h6th	 	 =		abs((r6th-R6th)*1.e-3)
h7th	 	 =		abs((r7th-R7th)*1.e-3)
h	 		 =	   np.concatenate((h1st ,h2nd, h3rd, h4th, h5th, h6th, h7th),axis =None)
h_asc		 =    np.concatenate((h1st, h2nd, h3rd, h4th, h5th, h6th),axis =None)

mt.plot(h)
mt.show()
#exit(0)
# RANGE ESTIMATION FROM VINCENTY FORMULA

Ran1st	 =		distVincenty(delta1st[0],delta1st[-1],lamda1st[0],lamda1st[-1])
Ran2nd	 =		distVincenty(delta1st[0],delta2nd[-1],lamda1st[0],lamda2nd[-1])
Ran3rd	 =		distVincenty(delta1st[0],delta3rd[-1],lamda1st[0],lamda3rd[-1])
Ran4th	 =		distVincenty(delta1st[0],delta4th[-1],lamda1st[0],lamda4th[-1])
Ran5th	 =		distVincenty(delta1st[0],delta5th[-1],lamda1st[0],lamda5th[-1])
Ran6th	 =		distVincenty(delta1st[0],delta6th[-1],lamda1st[0],lamda6th[-1])
Ran7th	 =		distVincenty(delta1st[0],delta7th[-1],lamda1st[0],lamda7th[-1])
RanApo	 =		distVincenty(delta1st[0],delta7th[0],lamda1st[0],lamda7th[0])

h_100 = np.where(h>100)[0]	# Flight in space, over 100 kilometers

over100_0 = h_100[0]
over100_f = h_100[-1]


microg_time = t[over100_f] - t[over100_0]

Q = dynamic_pressure(v,h) # Call the dynamic pressure function which gives the dynamic pressure [kPa]


Q_asc = dynamic_pressure(v_asc,h_asc)
Q_des = dynamic_pressure(v_des,h7th)



print('\n --------------FIRST FLIGHT PHASE-----------------------')
print('\n                     Powered                            ')
print(' Final speed at Talos burnout        = %f Km/s'% (v1st[-1]*1.e-3))
print(' Altitude at Talos burnout           = %f Km'% ( h1st[-1]))
print(' Range at Talos burnout              = %f Km'% ( Ran1st*1.e-3))
print(' Angle at Talos burnout              = %f deg'% ( phi1st[-1]/deg))
print(' Mass at Talos burnout               = %f Kg'% ( m1st[-1]))
print(' Propellant mass burned              = %f Kg'% ( mprop_talos))
print(' Talos ignition time                 = %f s'% ( tburn_talos))
print(' Mission elapsed time                = %f s'% ( tspan1[-1]))
print(' Drag loss				               = %f Km/s'% ( vlossdrag1st[-1]*1.e-3))
print(' Gravity loss				            = %f Km/s'% ( vlossgrav1st[-1]*1.e-3))
print('\n')
print('\n --------------SECOND FLIGHT PHASE----------------------')
print('\n                      Ballistic                         ')
print(' Mass at 1st Ballistic flight phase  = %f Kg'% ( m2nd[0]))
print(' Final speed at 1st Ballistic phase  = %f Km/s'% ( v2nd[-1]*1.e-3))
print(' Final altitude at 1st Ballis phase  = %f Km'% ( h2nd[-1]))
print(' Range after 1st Ballistic phase     = %f Km'% ( Ran2nd*1.e-3))
print(' Angle after 1st Ballistic phase     = %f deg'% ( phi2nd[-1]/deg))
print(' 1st Ballistic flight time           = %f s'% ( tspan2[-1]-tspan2[0]))
print(' Mission elapsed time                = %f s'% ( tspan2[-1]))
print(' Drag loss				               = %f Km/s'% ((vlossdrag2nd[-1]-vlossdrag2nd[0])*1e-3))
print(' Gravity loss				            = %f Km/s'% ( (vlossgrav2nd[-1]-vlossgrav2nd[0])*1e-3))
print('\n')
print('\n --------------THIRD FLIGHT PHASE-----------------------')
print('\n                     Powered                            ')
print(' Final speed at Tomahawk burnout     = %f Km/s'% ( v3rd[-1]*1.e-3))
print(' Altitude at Tomahawk burnout        = %f Km'% ( h3rd[-1]))
print(' Range at Tomahawk burnout           = %f Km'% ( Ran3rd*1.e-3))
print(' Angle at Tomahawk burnout           = %f deg'% ( phi3rd[-1]/deg))
print(' Mass at Tomahawk burnout            = %f Kg'% ( m3rd[-1]))
print(' Propellant mass burned              = %f Kg'% ( mprop_toma))
print(' Tomahawk ignition time              = %f s'% ( tburn_toma))
print(' Mission elapsed time                = %f s'% ( tspan3[-1]))
print(' Drag loss                        	   = %f Km/s'% ((vlossdrag3rd[-1]-vlossdrag3rd[0])*1.e-3))
print(' Gravity loss				            = %f Km/s'% ( (vlossgrav3rd[-1]-vlossgrav3rd[0])*1.e-3))
print('\n')
print('\n --------------FOURTH FLIGHT PHASE----------------------')
print('\n                     Ballistic                          ')
print(' Mass at 2nd Ballistic flight phase  = %f Kg'% ( m4th[0]))
print(' Final speed at 2nd Ballistic phase  = %f Km/s'% ( v4th[-1]*1.e-3))
print(' Final altitude at 2nd Ballis phase  = %f Km'% ( h4th[-1]))
print(' Range after 2nd Ballistic phase     = %f Km'% ( Ran4th*1.e-3))
print(' Angle after 2nd Ballistic phase     = %f deg'% ( phi4th[-1]/deg))
print(' 2nd Ballistic flight time           = %f s'% ( tspan4[-1]-tspan4[0]))
print(' Mission elapsed time                = %f s'% ( tspan4[-1]))
print(' Drag loss				               = %f Km/s'% ( (vlossdrag4th[-1]-vlossdrag4th[0])*1.e-3))
print(' Gravity loss				            = %f Km/s'% ( (vlossgrav4th[-1]-vlossgrav4th[0])*1.e-3))
print('\n')
print('\n --------------FIFTH FLIGHT PHASE-----------------------')
print('\n                     Powered                            ')
print(' Final speed at Nikha burnout        = %f Km/s'% ( v5th[-1]*1.e-3))
print(' Altitude at Nikha burnout           = %f Km'% ( h5th[-1]))
print(' Range at Nikha burnout              = %f Km'% ( Ran5th*1.e-3))
print(' Angle at Nikha burnout              = %f deg'% ( phi5th[-1]/deg))
print(' Mass at Nikha burnout               = %f Kg'% ( m5th[-1]))
print(' Propellant mass burned              = %f Kg'% ( mprop_nikha))
print(' Nikha ignition time                 = %f s'% ( tburn_nikha))
print(' Mission elapsed time                = %f s'% ( tspan5[-1]))
print(' Drag loss				               = %f Km/s'% ( (vlossdrag5th[-1]-vlossdrag5th[0])*1.e-3))
print(' Gravity loss				            = %f Km/s'% ( (vlossgrav5th[-1]-vlossgrav5th[0])*1.e-3))
print('\n')
print('\n --------------SIXTH FLIGHT PHASE-----------------------')
print('\n                     Ballistic                          ')
print(' Mass at 3rd Ballistic flight phase  = %f Kg'% ( m6th[0]))
print(' Final speed at 3rd Ballistic phase  = %f Km/s'% ( v6th[-1]*1.e-3))
print(' Final altitude at 3rd Ballis phase  = %f Km'% ( h6th[-1]))
print(' Range after 3rd Ballistic phase     = %f Km'% ( Ran6th*1.e-3))
print(' Angle after 3rd Ballistic phase     = %f deg'% ( phi6th[-1]/deg))
print(' 3rd Ballistic flight time           = %f s'% ( tspan6[-1]-tspan6[0]))
print(' Mission elapsed time                = %f s'% ( tspan6[-1]))
print(' Drag loss				               = %f Km/s'% ( (vlossdrag6th[-1]-vlossdrag6th[0])*1.e-3))
print(' Gravity loss				            = %f Km/s'% ( (vlossgrav6th[-1]-vlossgrav6th[0])*1.e-3))
print('\n')
print('\n ---------------SEVENTH FLIGHT PHASE--------------------')
print('\n                     Reentry                            ')
print(' Final speed                         = %f Km/s'% ( v[-1]))
print(' Speed at apogee                     = %f Km/s'% ( v7th[0]*1.e-3))
print(' Apogee Altitude                     = %f Km'% ( max(h)))
print(' Apogee Range                        = %f Km'% ( RanApo*1.e-3))
print(' Total Range                         = %f Km'% ( Ran7th*1.e-3))
print(' Final angle                         = %f deg'% ( phi[-1]))
print(' Final Mass                          = %f Kg'% ( m[-1]))
print(' Reentry flight time                 = %f s'% ( tspan7[-1]-tspan7[0]))
#print(' Microgravity time                   = %f s'% ( microg_time))
print(' Total flight time                   = %f s'% ( t[-1]))
print(' Final altitude                      = %f Km'% ( h[-1]))
print(' Total Drag loss (only ascent)             = %f Km/s'% ( vlossdrag[-1]))
print(' Total Gravity loss (only ascent)          = %f Km/s'% ( vlossgrav[-1]))
print('\n\n ---------------------------------------------------\n')





# figure(1) # Plots the altitude as a function of time
# plot(t,h)
# axis equal
# xlabel('Time (s)')
# ylabel('Altitude (km)')
# axis([-inf, inf, -inf, inf])
# grid
# 
# 
# figure(2) # Plots the velocity as a function of time
# plot(t, v)
# xlabel('time (s)')
# ylabel('velocity (km/s)')
# grid
# 
# figure(3) # Plots the angle wrt local horizon as a function of time
# plot(t,phi)
# xlabel ('Time')
# ylabel ('Angle (deg)')
# grid
# 
# figure(4) # Plots the mass as a function of time
# plot(t,m)
# xlabel('Time (s)')
# ylabel('Mass (Kg)')
# axis([0, 100, 0, inf])
# grid
# 
# figure(5) # Plots the Ground Track
# plot(lamda,delta,'g')
# xlabel('Longitude')
# ylabel('Latitude')
# 
# figure(6) # Plots the velocity as a function of altitude
# plot(h_asc,v_asc,'r',h7th,v_des,'b')
# xlabel('altitude (km)')
# ylabel('velocity (km/s)')
# grid
# 
# figure(7) # Plots the Dynamic pressure as a function of altitude
# plot(Q_asc,h_asc,'r',Q_des,h7th,'g')
# xlabel('Dynamic pressure (kPa)')
# ylabel('Altitude (Km)')
# axis([-inf, inf, -inf, 80])
# grid
# 
# figure(8) # Plots the dynamic pressure as a function of time
# plot(t,Q,'g')
# xlabel('Time')
# ylabel('Dynamic pressure ascent (kPa)')
# axis([-inf, 60, -inf, inf])
# grid





tiempo_over100 = t[over100_0:over100_f+1]
delta_over100= delta[over100_0:over100_f+1]*deg
lamda_over100=lamda[over100_0:over100_f+1]*deg
phi_over100=phi[over100_0:over100_f+1]*deg
v_over100=v[over100_0:over100_f+1]
Azimut_over100=A[over100_0:over100_f+1]*deg
radio_over100 = r[over100_0:over100_f+1]*1.e+3
h_over100 = h[over100_0:over100_f+1]

#x,y,z coordenadas trayectoria total

x_t = (r*1e+3)*np.cos(delta*deg)*np.cos(lamda*deg)
y_t = (r*1e+3)*np.cos(delta*deg)*np.sin(lamda*deg)
z_t = (r*1e+3)*np.sin(lamda*deg)

#coordenadas en metros 
x = radio_over100*np.cos(delta_over100)*np.cos(lamda_over100) # delta = latitud , lamda = longitud 
y = radio_over100*np.cos(delta_over100)*np.sin(lamda_over100)
z = radio_over100*np.sin(lamda_over100) 

# coordenadas en km 
x_2 = (radio_over100*np.cos(delta_over100)*np.cos(lamda_over100))*1.e-3 
y_2 =( radio_over100*np.cos(delta_over100)*np.sin(lamda_over100))*1.e-3
z_2 =( radio_over100*np.sin(delta_over100))*1.e-3 

#Velocidades VxNorth,VyEast, VzDown
vNorth = v_over100*np.cos(Azimut_over100)*np.cos(phi_over100)
vEast = v_over100*np.cos(phi_over100)*np.sin(Azimut_over100)
vDown= -1*v_over100*np.sin(phi_over100)

#Vx,Vy,Vz
[Vx,Vy,Vz] = ned2ecefv(vNorth, vEast, vDown, delta_over100, lamda_over100)
#[Vx,Vy,Vz]=ned2ecefv(vNorth,vEast,vDown,delta_over100,lamda_over100,'radians')


Tabelle ={'Posicion':h_100,
          'Tiempo[s]':tiempo_over100,
          'Latitud[rad]':delta_over100,
          'Longitud[rad]':lamda_over100,
          'Angulo_horizonte[rad]':phi_over100,
          'Velocidad[km/s]':v_over100,
          'Azimut[rad]':Azimut_over100,
          'Altura[km]':h_over100,
          'Radio[m]':radio_over100,
          'X':x_2,
          'Y':y_2,
          'Z':z_2,
          'North[km/s]':vNorth,
          'East[km/s]':vEast,
          'Down[km/s]':vDown,
          'Vx[km/s]':Vx,
          'Vy[km/s]':Vy,
          'Vz[km/s]':Vz
          }
df = DataFrame(Tabelle, columns= ['Posicion','Tiempo[s]','Latitud[rad]','Longitud[rad]',
                                  'Angulo_horizonte[rad]','Velocidad[km/s]',
                                  'Azimut[rad]','Altura[km]','Radio[m]','X','Y','Z',
                                  'North[km/s]','East[km/s]','Down[km/s]','Vx[km/s]',
                                  'Vy[km/s]','Vz[km/s]'])


export_excel = df.to_excel (r'./export_dataframe.xlsx',
                            index = None, header=True)         
          
          
          
          
          






