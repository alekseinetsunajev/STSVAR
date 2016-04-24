%Function to calculate likelihood
function [L logLikVec] = calculate_L(obj)

u = residuals2(obj.T, obj.Y, obj.Z, obj.lags, obj.Theta);

% calculate the initial likelihood value for each t.

for t = 1: obj.T(1,1)-obj.lags    
       Exponent = u(t,:)  * (((1- obj.TR_f(t,1))*(obj.B*obj.B') + obj.TR_f(t,1) * (obj.B*obj.Lambda*obj.B'))^-1 )* u(t,:)';   
       Det = (det((1-obj.TR_f(t,1))*(obj.B*obj.B') + obj.TR_f(t,1) * (obj.B*obj.Lambda*obj.B'))) ;
       logLikVec(t,1) = log(sqrt(1/((2*pi)^obj.T(1,2)))*1/sqrt(Det) * exp(-0.5 * Exponent));
end
            
L = sum(logLikVec);