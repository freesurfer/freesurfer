function h = fast_hann(x,tau)
% h = hann(x,tau)
%
% Hann window 
%   alpha = 0.5; % 0.54 for Hamming
%   h = alpha + (1-alpha)*cos(pi*x(ind)/tau);
%
% $Id: fast_hann.m,v 1.1 2007/01/14 05:39:07 greve Exp $

alpha = 0.5; % 0.54 for Hamming
h = zeros(size(x));
ind = find(abs(x) < tau);
h(ind) = alpha + (1-alpha)*cos(pi*x(ind)/tau);

return;
