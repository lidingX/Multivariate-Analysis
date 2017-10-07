function model = logReg(X, Y, epsilon,  t)
% Input:
%   X: n x d data maxtrix
%   Y: n x 1 label (0/1)
%   epsilon: convergence tolerence
%   t : number of iteration
% Output:
%   model: model of logistic regression
% Written by LiDing Xu(a648219931@163.com).    
if(nargin < 4)
    t = 1e5; 
elseif(nargin < 3)
    epsilon = 1e-4;
end
[n,~] = size(Y);
X = [ones(n,1) X];
z = log((Y+0.5)./(1.5-Y));   
v = 6./((Y+1).*(2-Y));
D = diag(1./v);
beta = (X'*D*X)\X'*D*z; %Modification of the Empirical Logistic Transformation
N = exp(X*beta)./(1+exp(X*beta));
A = diag(N.*(1-N));
q = X'*(Y - N);
H = -X'*A*X;
while(norm(q)>epsilon && t < 1e6) %Newton-Raphson iteration
    beta = beta - H\q;
    N = exp(X*beta)./(1+exp(X*beta));
    A = diag(N.*(1-N));
    q = X'*(Y - N);
    H = -X'*A*X;
    t = t + 1;
end
beta = beta - H\q;
model.beta = beta;
end

function P = logRegPred(X,model)
% Input:
%   X: n x d data maxtrix
%   model: model of logistic regression
% Output:
%   P: probability of logistic regression
% Written by LiDing Xu(a648219931@163.com).  
P = exp(X*model.beta)./(1+exp(X*model.beta));
end