function [xp] = mathematicalmodelukf(x)
xp(1) = x(4)*(x(2)-x(1));
xp(2) = x(1)*(x(5)-x(3))-x(2);
xp(3) = x(1)*x(2)-x(6)*x(3);
xp(4) = 0;
xp(5) = 0;
xp(6) = 0;
end