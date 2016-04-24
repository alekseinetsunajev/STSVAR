function [LogL] = minimize(Y, Z, T, lags, Param, Theta, TR_f)

u = residuals2(T, Y, Z, lags, Theta);

B = reshape(Param( 1 :   T(1,2)^2), T(1,2),  T(1,2));
L = diag(Param( 1+ T(1,2)^2 : T(1,2)^2 + T(1,2)));

logDet = 0;
for t = 1: T(1,1) - lags    
    R(t, 1) = u(t,:)  * (((1 - TR_f(t,1))*(B*B') + TR_f(t,1) *(B*L*B'))^-1 )* u(t,:)';   
     %(1-exp(-1*gamma * (t - c)^2))
     logDet = logDet + log(det((1- TR_f(t,1))*(B*B') + TR_f(t,1)* (B*L*B'))) ;
end
SSR = sum(R);         

if isreal( 0.5*logDet + 0.5 *SSR)
    LogL = 0.5*logDet + 0.5 *SSR;
else
    B;
    SSR;
    LogL = 1000;
end