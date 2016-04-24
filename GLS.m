function [Theta] = GLS(Z, Y, B, L, T, lags, TR_f)

for t = 1:T(1,1) - lags
      Omega((t-1)*T(1,2) + 1 : t*T(1,2), (t-1)*T(1,2) +1 : t*T(1,2)) = ( (1 - TR_f(t,1))*(B*B') + ...
                    TR_f(t,1) * (B*L*B'))^-1;
end

V = (kron(Z', eye(T(1,2))) * Omega *  kron(Z, eye(T(1,2))))^-1 ;
y = reshape(Y', T(1,2)*T(1,1), 1);

% update parameters
Theta = V * kron(Z', eye(T(1,2))) * Omega * y(T(1,2)*lags+1:end, :);