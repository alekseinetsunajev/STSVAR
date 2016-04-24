% function to obtain the impulse responses
function [IR_SR, IR_LR]=irs(obj)

% rescale B such that it has ones on the main diagonal
A = zeros(obj.T(1, 2)* (obj.lags) ); % Define the A matrix
for i = 1:obj.lags
       A(1:obj.T(1,2), obj.T(1,2)*(i-1)+1 : obj.T(1,2)*i ) = get_coefficient(obj.Theta, obj.T, i, obj.trend);
end
for m=1:obj.lags-1
       A(m*obj.T(1,2)+1 : m*obj.T(1,2)+obj.T(1,2), obj.T(1,2)*m-obj.T(1,2)+1 : obj.T(1,2) * m) = eye(obj.T(1,2));
end
J = [eye(obj.T(1,2)) zeros(obj.T(1,2), obj.T(1,2)*obj.lags - obj.T(1,2) )  ]; % Define the J matrix

% short run responses
IR_SR = reshape(J*A^0*J'*obj.B, obj.T(1,2)^2, 1); % Impulse response matrix
for i=1:obj.h
    IR_SR=([IR_SR reshape(J*A^i*J'*obj.B, obj.T(1,2)^2, 1)]);
end;
            
% long run responses
IR_LR = IR_SR(:,1);
for i = 1 : obj.h
      IR_LR(:,i+1) =  IR_LR(:, i) + IR_SR(:, i+1);
end
            
