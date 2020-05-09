function [tmean,tsigmapoints,tP,tdeviations]=uth(X,Wm,Wc,R)
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
L     = size(X,1);
tmean = zeros(1,2);

Y     = X(:,1:2);
tmean = tmean+Wm*Y;

tsigmapoints = Y;
tdeviations  = Y-tmean(ones(1,L),:);
tP=tdeviations'*diag(Wc)*tdeviations+R;
end