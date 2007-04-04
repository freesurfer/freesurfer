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
% dlh = 4.6167/(prod(fwhm))
% 
% $Id: fast_dlh.m,v 1.2 2007/04/04 02:07:46 greve Exp $

gstd = fast_fwhm2std(fwhm);
dlh = 1/(sqrt(8)*prod(gstd));

return;





