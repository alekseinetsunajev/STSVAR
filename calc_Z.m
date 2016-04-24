function [Z] = calc_Z (y, lags)

T = size(y);
Z = zeros(T(1,1)- lags, lags *T(1,2) + 1 );
for i = 1:T(1,1)- lags;
    Z(i,1) = 1;                                         % constant
    for j=1: lags
        Z(i,(j-1)*T(1,2)+2 : j*T(1,2)+1 )= y(lags+i-j,:);  % lagged observations 
    end
end