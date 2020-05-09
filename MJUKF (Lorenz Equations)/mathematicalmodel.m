function [xp] = mathematicalmodel(x, theta)
global time
% Step change in real parameters
% if time>40
%     theta(1)=7;
%     theta(2)=21;
%     theta(3)=1;
% end
xp(1) = theta(1)*(x(2)-x(1));
xp(2) = x(1)*(theta(2)-x(3))-x(2);
xp(3) = x(1)*x(2)-theta(3)*x(3);
end