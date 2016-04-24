function[C] = get_coefficient(Theta, T, Position, trend)

if trend ==1
     cnt1 = (Position-1)*T(1,2)^2 + 2*T(1,2) + 1;   
     cnt2 = (Position)*T(1,2)^2 + 2*T(1,2);  
else
    cnt1 = (Position-1)*T(1,2)^2 + T(1,2) + 1;   
    cnt2 = (Position)*T(1,2)^2 + T(1,2);
end
C = reshape(Theta(cnt1:cnt2), T(1,2), T(1,2) );