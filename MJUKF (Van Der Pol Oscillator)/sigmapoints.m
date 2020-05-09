function X=sigmapoints(x,P,gamma)
%Sigma points around reference point
%Inputs: 
%       x:     reference point
%       P:     covariance
%       gamma: scaling parameter
%Output:
%       X: Sigma points
A = gamma*chol(P)';
Y = x(ones(1,numel(x)),:);
X = [x; Y+A; Y-A]; 
end