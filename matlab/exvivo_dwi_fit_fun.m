function yhat = exvivo_dti_fit_fun(coeffs, x)
% coeffs are the coefficients to fit
% x is a vector of flip angles
% y (yhat) is the found (estimated) signal intensity for a given flip angle
S0 = coeffs(1);
tau = coeffs(2);
delta = coeffs(3);
offset = coeffs(4);

yhat = S0 * exp(-(x-delta)/tau) + offset;
