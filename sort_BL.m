function[B L] = sort_BL(Param, T)
% sort Lambda in increasing order and adjust B accordingly
% Param  - [vec(B); diag(Lambda)]
B = reshape(Param( 1 :   T(1,2)^2), T(1,2),  T(1,2));
L = diag(Param( 1+ T(1,2)^2 : T(1,2)^2 + T(1,2)));

L0 = length(L);
for i = 1: L0-1
        for j=i+1:L0
            if L(i,i) > L(j,j)
                L_temp = L(i,i);
                L(i,i) = L(j,j);
                L(j,j) = L_temp;

                B_temp(:,1) =  B(:,i);
                B(:,i) = B(:,j);
                B(:,j) = B_temp(:, 1);
            end
        end
end
    
 %        Check signs of B
            if B(1,1) < 0
                B(:,1) = -1 * B(:,1);
            end
            if B(2,2) < 0
                B(:,2) = -1 * B(:,2);
            end
            if B(4,3) < 0
                B(:,3) = -1 * B(:,3);
            end
            if B(3,4) < 0
                B(:,4) = -1 * B(:,4);
            end
            if B(5,5) < 0
                B(:,5) = -1 * B(:,5);
            end            
