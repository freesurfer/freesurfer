function [eres, sigest] = fmri_residual(y,X,hest)
%
% [eres sigest] = fmri_residual(y,X,hest)
%
% Computes the residual error between the actual signal y
% and the signal estimate sigest = X*h.
%
% $Id: fmri_residual.m,v 1.1 2003/03/04 20:47:40 greve Exp $

sigest = fmri_estsignal(X,hest);
eres = y - sigest;

return;
