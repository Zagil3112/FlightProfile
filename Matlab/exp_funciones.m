 
%Funciones independientes 

x_exp4 = atmosphere(0.01,2000);
Cd_exp = Cd_reentry1221(0.001);  

%Funciones independientes 
x_exp1 = powered(0,z1_0,p1);
% x_exp2 = ballistic(0,z2_0);
% x_exp3 = reentry(0,z7_0);

%Integrador 
% [t1,z1]=ode45(@powered,tspan1,z1_0,[],p1);
% [t2,z2]=ode45(@ballistic,tspan2,z2_0);
[t7,z7]=ode45(@reentry,tspan7,z7_0,options2);


% dist_exp = distVincenty(42.3541165,40.7791472,-71.0693514,-73.9680804);
% a_zap = [2;4 ;6];
% b_zap = [ 3 ;5 ;7];
% c_zap = [ 9 ;9 ;9];
% 
% d_zap = [a_zap; b_zap; c_zap]
% 
% v_exp = v(1:6)
% h_exp = h(1:6)







    












