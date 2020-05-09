function [xcap,P,thetacap,thetacapvec] = unscentedkalmanfilter(xcap,P,x,thetacap,thetacapvec)
global gamma
global Wm
global Wc
global nse
global Q
global R

%sigma points around x
X = sigmapoints(xcap,P,gamma); 
[txcap,tX,transformedP,Xdev,~] = utf(X,Wm,Wc,nse,Q,thetacapvec); 
[txcap1,tY,transformedP1,Xdev1] = uth(tX,Wm,Wc,nse,R); 

P12  = Xdev'*diag(Wc)*Xdev1;   %transformed cross-covariance
K    = P12*inv(transformedP1);
xcap = txcap+(K*(x-txcap1'))'; %state update
P    = transformedP-K*P12';    %covariance update

thetacapvec1 = thetacap - 0.25*[0.5,0.5]*(x-tY');
thetacapvec  = thetacapvec1';
thetacap     = mean(thetacapvec);
end