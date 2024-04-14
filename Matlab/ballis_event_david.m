function[lookfor stop direction]=ballis_event(t,x)

global Re epsilon;

lookfor=x(1); %searches for this expression set to 0
stop=1; %stop when event is located.
direction=-1; %specify direction of motion at event
