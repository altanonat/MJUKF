function [xp] = mathematicalmodelukf(x,theta)
xp(1) = x(2);
xp(2) = theta*(1-(x(1))^2)*x(2)-x(1);
end