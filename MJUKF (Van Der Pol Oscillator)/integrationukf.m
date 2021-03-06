function [x] = integrationukf(x,xp,theta)
global deltat
global time 

tempx = x;
k1    = deltat*xp;
x     = tempx+0.5*k1;
time  = time+0.5*deltat;
xp    = mathematicalmodelukf(x,theta);

k2 = deltat*xp;
x  = tempx+0.5*k2;
xp = mathematicalmodelukf(x,theta);

k3   = deltat*xp;
x    = tempx+k3;
time = time+0.5*deltat;
xp   = mathematicalmodelukf(x,theta);

k4 = deltat*xp;
x  = tempx+(1/6)*(k1+2*k2+2*k3+k4);
end