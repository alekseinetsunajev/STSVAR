% This program create a vector from a matrix
% V = vec(X,rc)
% X must be a matrix 
% stacking rows rc=1 or columns rc=2.
function V = vec(X,rc)
[a b] = size(X);
if rc==1
    for i=1:a
        V(1,(i-1)*b+1:i*b) = X(i,:);  %% stack rows
    %V((i-1)*a+1:i*a,1)=X(:,i);      %% stack columns
    end
    V = V';
elseif rc==2
    for i=1:b
        V((i-1)*a+1:i*a,1)=X(:,i); %% stack columns
    end
end
