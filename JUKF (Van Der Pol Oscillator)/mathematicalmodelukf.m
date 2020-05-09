function [xp] = mathematicalmodelukf(x)
xp(1) = x(2);
xp(2) = x(3)*(1-(x(1))^2)*x(2)-x(1);
xp(3) = 0;
end