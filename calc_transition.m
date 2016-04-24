function tr_value = calc_transition(t, c, gamma, exponent)

if exponent == 1
    tr_value = (1+exp(-1*exp(gamma) * (t - c)))^-1;
else
    tr_value = (1+exp(-1*gamma * (t - c)))^-1;
end