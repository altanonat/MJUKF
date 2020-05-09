global alpha
global kappa
global beta
global gamma
global lambda
global Px
global q
global r
global nse
global Wm
global Wc
global Q
global R

alpha = 1;    %default, tunable
kappa = 0;    %default, tunable
beta  = 2;    %default, tunable

lambda = alpha^2*(nse+kappa)-nse; %scaling factor
gamma  = nse+lambda; %scaling factor

%weights for means
Wm    = [lambda/gamma (0.5)./(gamma+zeros(1,2*nse))];
Wc    = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta); %weights for covariance
gamma = sqrt(gamma);

q  = 1e-3;
r  = 1e-1;

Px = 5*eye(nse);

Q  = q*eye(nse);
R  = r*eye(nse);