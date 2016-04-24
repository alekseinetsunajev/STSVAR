function VarDec = FEVD(obj, State)
%horizon = 100;
% calculate positions of IRs, example {1 4 7 2 5 8 . 6 9 }
Vec = ones(obj.T(1,2)^2, 1 );
for i = 2: obj.T(1,2)^2
    Vec(i) = Vec(i-1) + 1;
end
Vec = reshape ( reshape(Vec, obj.T(1,2), obj.T(1,2))', obj.T(1,2)^2, 1 )';          

if State == 1
    L = ones(obj.T(1,2), 1);
    L_dup = kron(sqrt(L), ones(obj.T(1,2), 1)); 
else
    L = diag(obj.Lambda);                                                               % diagonal of Lambda 
    L_dup = kron(sqrt(L), ones(obj.T(1,2), 1));                              % weight of the Lambda
end
% MSE computation
MSE=zeros(obj.T(1,2) , obj.h);
for j = 1:obj.h
    for i = 1: obj.T(1,2)   % var
       Mse = 0;
       for k = 1:obj.T(1,2)
            Mse = Mse +  ( ( sqrt(L(k)) * obj.IR_SRun( Vec( (i-1)*obj.T(1,2) + k ), j) ) ^2 ) ;
       end
       if j > 1
            MSE(i, j) = MSE(i, j-1) +  Mse;
       else
           MSE(i, j) =  Mse;
       end     
    end
end

% computation of each component
Theta_Sq(:, 1) = ( (L_dup .* obj.IR_SRun(:, 1) ) .^ 2 ) ;
for j = 2:obj.h
    for i = 1: obj.T(1,2)^2   % var
            Theta_Sq(i, j) =  Theta_Sq(i, j-1) + ( (L_dup(i) * obj.IR_SRun(i, j) )  .^2)   ;
     end
end

% Computation of the contribution of each shock to total variance
for  j = 1:obj.h
    cnt = 1;
    for i = 0 : obj.T(1, 2) : ( obj.T(1, 2)^2 - 1 ); 
        for k = 1:obj.T(1,2)
            VarDec(i+k, j) = Theta_Sq( Vec(i+k), j) / MSE(cnt, j);
        end
        cnt = cnt + 1;
    end
end



