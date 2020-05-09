global alpha
global kappa
global beta
global gamma
global lambda
global Px
global q
global r
global qp
global sigma
global rho
global zheta
global nse
global Wm
global Wc
global Q
global R

alpha = 1;    %default, tunable
kappa = 0;    %default, tunable
beta  = 2;    %default, tunable

sigma = 10;
rho   = 28;
zheta = 8/3;

lambda = alpha^2*(nse+kappa)-nse; %scaling factor
gamma  = nse+lambda; %scaling factor

%weights for means
Wm    = [lambda/gamma (0.5)./(gamma+zeros(1,2*nse))];
Wc    = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta); %weights for covariance
gamma = sqrt(gamma);

q  = 1e-3;
r  = 1e-1;
qp = 0.5;

Px          = q*eye(nse);
Px(4:6,4:6) = qp*eye(3);

Q  = q*eye(nse);
R  = r*eye(3);