function [tmean,tsigmapoints,tP,tdeviations]=utf(X,Wm,Wc,n,R,yms)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed sampling points
%        P: transformed covariance
%       Y1: transformed deviations
global deltat
global time 
L     = size(X,1);
tmean = zeros(1,n);
Xdot  = zeros(L,n);
for k=1:L 
    Xdot(k,:) = mathematicalmodelukf(X(k,:));
    X(k,:)    = integrationukf(X(k,:),Xdot(k,:)); 
    time      = time-deltat;
    tmean     = tmean+Wm(k)*X(k,:);       
end
tsigmapoints = X;
tdeviations  = X-tmean(ones(1,L),:);
tP=tdeviations'*diag(Wc)*tdeviations+R;
end