function acf = fast_fbmacf(H,nmax)
% acf = fast_fbmacf(H,nmax)
% 
% Autocorrelation function for fractional Brownian motion
%
% Fadili and Bullmore, Wavelet-Generalized Least Squares,
% NI 15, 217-232, 2002. Equation 2.
%
% H = 0.5, process is white
% .5 < H < 1, process has long-range dependence
% 0 < H < .5, process has short-range dependence and
%  will "anti-persistence", negative correlation
%
% nmax - max number of lags
%

lag = [0:nmax-1]';

acf = ((lag+1).^(2*H) - 2*lag.^(2*H) + (lag-1).^(2*H))/2;
ind0 = find(lag == 0);
acf(ind0) = 1;

return;





