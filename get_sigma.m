% Retrieve the covariance matrix for prespecified position
function[Sigma] = get_sigma(input_sigma, T, position)

Sigma = input_sigma (1 + T(1,2)*(position-1) : T(1,2) * position, :) ;

