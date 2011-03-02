function [havg, hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% [havg hstd] = fmri_avgrandeff(hsum, hsum2, dof)
%
% Comptutes the average and std dev for rand effects accumulated
% values.
%
%


%
% fmri_avgrandeff.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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
