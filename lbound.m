% function returns the matrix of lower bounds for the estimates of B and
% Lambda
function[lb]=lbound(T, s)

for a = 1:T(1,2)
    for b = 1:T(1,2)
        lb(a,b) = -100;
    end
end

for a = 1:s-1
    lb(:, T(1,2)+a) = 0.000001*ones(1, T(1,2))';
end
