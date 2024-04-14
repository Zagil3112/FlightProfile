function[lookfor stop direction]=reentry_event(t,x)

global Re epsilon;

R=Re*(1-epsilon*(sin(x(5)))^2)+5000;

lookfor=x(4)-R; %searches for this expression set to 0
stop=1; %stop when event is located.
direction=-1; %specify direction of motion at event
