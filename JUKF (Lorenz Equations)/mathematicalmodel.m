function [xp] = mathematicalmodel(x)
global sigma
global rho
global zheta
global time

% if time>10
%     sigma = 7;
%     rho   = 21;
%     zheta = 1;
% end

xp(1) = sigma*(x(2)-x(1));
xp(2) = x(1)*(rho-x(3))-x(2);
xp(3) = x(1)*x(2)-zheta*x(3);
end