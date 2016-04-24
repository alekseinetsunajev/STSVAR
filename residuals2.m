function[u] = residuals2(T, y, Z, lags, Theta)
u = reshape(reshape(flipud(rot90(y(lags+1:T(1,1), :) )), (T(1,1)-lags)*T(1,2), 1) - kron(Z, eye(T(1,2)))*Theta, ...
    T(1,2), T(1,1)-lags)';