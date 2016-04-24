function [LogL] = LogLikeVAR(vB, T, lags, U)
B = reshape(vB, T(1,2), T(1,2));        
SSR = trace(U * (B*B')^-1 * U');
logDet = (T(1,1) - lags) * log(det(B*B' ));

if isreal( 0.5*logDet + 0.5 *SSR)
    LogL =   0.5*logDet + 0.5 *SSR;
else
    B
    SSR
    LogL = 1000;
end
