function [havg, hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% [havg hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% Comptutes the average and std dev for rand effects accumulated
% values.
%
% $Id: fmri_avgrandeff.m,v 1.1 2003/03/04 20:47:39 greve Exp $

if(nargin ~= 3)
  msg = 'USAGE: [havg hstd] = fmri_avgrandeff(hsum, hsum2, dof)'
  qoe(msg);error(msg);
end

if(dof < 2)
  msg = sprintf('dof = %d, must be > 1',dof);
  qoe(msg);error(msg);
end

havg = hsum/dof;
hstd = real(sqrt( (hsum2 - dof*(havg.^2))/(dof-1)));

%tmp = (hsum2 - dof*(havg.^2))/(dof-1);
%keyboard
return;