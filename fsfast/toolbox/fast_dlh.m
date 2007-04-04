function dlh = fast_dlh(fwhm)
% dlh = fast_dlh(fwhm)
%
% fwhm is the fwhm for each dim
%
% 1/sqrt(det(Lambda))
% This is what FSL uses
%
% gstd = fast_fwhm2std(fwhm);
% dlh = 1/(sqrt(8)*prod(gstd));
% 
% $Id: fast_dlh.m,v 1.1 2007/04/04 02:00:08 greve Exp $

gstd = fast_fwhm2std(fwhm);
dlh = 1/(sqrt(8)*prod(gstd));

return;





