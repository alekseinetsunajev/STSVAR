%% OLS ESTIMATES OF A VAR(p) MODEL, conditioned on the first p observations
% MODEL: y_t=v+A1*y_t-1+...+Ap*y_t-p+u_t, 
% where y_t is a kx1 vector, each Ai a KxK matrix (i=1,...,p) and u_t a Kx1
% vector.
% INPUTS: %%
% y - a KxT matrix of data, K num. of var., T num. of obs.
% p - lag length of VAR model
% OUTPUTS: %%
% K - number of eq.
% T - number of obs.
% p - selected lag length
% beta - a Kx(1+Kp) matrix of estimated coeff., beta=[v A1 ... Ap], where v
% is a Kx1 vector of constants, Ai a KxK matrix of i-th lag estimated coeff
% and i=1,...,p
% cvm - a KxK var-cov. matrix of est. resid.
% sbeta - a Kx(1+Kp) matrix of std. err. of the est. coeff.
% tbeta - a Kx(1+Kp) matrix of t stat. of the est. coeff
% resids - est. residuals

function [beta, resids, cvm, SSR]=Var(y,p)


%% PRELIMINARIES %%
% get dimensions
[K T]=size(y);

%% TRANSFORMATION of THE MODEL %%
% Into the following form: Y=beta*Z+U, where:
% Y:=[y_p+1 ... y_T] is a Kx(T-p) matrix, we condition on the first p obs.
% beta:= [v A1 ... Ap], as before
% Z:= [Z_p ... Z_T-1] is a (1+pK)x(T-p+1) matrix, where Z_t=[1 y_t ...
% y_t-p+1]', y_t as before, t=p,...,T-1 (see Lutkepohl, 2007, p.70)
% U:=[u_p+1 ... u_T] is a Kx(T-p) matrix of est. resid.


%% Builds the Y matrix
Y=y(:,p+1:T);

%% Builds the Z matrix
Z=zeros(1+K*p,T-p);
Z(1,:)=1;
for i=0:p-1
    Z(i*K+2:(i+1)*K+1,:)=y(:,p-i:T-i-1);
end

%% ESTIMATES %%
% creates beta mat. (see Lutkepohl, 2007, p.72, formula 3.2.10)
B=(Y*Z')/(Z*Z'); 
beta = reshape(B, K + p*K^2, 1);
%% creates residuals
resids = Y - B*Z;
% creates covar mat. of resid. (See Lutkepohl, 2007, p.75, formula 3.2.18 and 3.2.19)
cvm=resids*resids'*1/(K*(T - p) - (K*p+1)*K);

cov_beta = kron((Z*Z')^-1, cvm);
SSR = resids*resids';
%% Creates sbeta and tbeta 
% creates a Kx(1+Kp) matrix sbeta 
% first we create a var-cov matrix of "vec(beta)" and call it cVecBeta (see Lutkepohl, 2007, p. 77, formula 3.2.21) 
%[estimates of std. err. of vec(beta) are then the corresponding square roots of the diagonal elements of cVecBeta]

%cVecBeta=kron(inv(Z*Z'),cvm); 
% calculates std. errors of vec(beta) by retriewing diagonal elements of
% cVecBeta and takes their sqrt.
%sVecBeta=sqrt(diag(cVecBeta)); 
% reverses the vectorized sVecBeta matrix into sbeta matrix of std. err. of 
% the corresponding est. coeff. 
%sbeta=reshape(sVecBeta,K,K*p+1); 
% creates tbeta
%tbeta=beta./sbeta; %CHECK!!

%% Clears the work space of unneccessary results
%clear sVecBeta sVecBeta Z cVecBeta i;

%% ADD THE FOLLOWING STATISTICS

%% R^2
%rsqr = 1.0 - SSR/(y - mean(y))'*(y - mean(y)) % Check leSage, eq. by eq
% rbar, adjusted 
% D-W stat
% Q stat
















