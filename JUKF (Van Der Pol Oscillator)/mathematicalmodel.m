function [xp] = mathematicalmodel(x,mu)
xp(1) = x(2);
xp(2) = mu*(1-(x(1))^2)*x(2)-x(1);
end
