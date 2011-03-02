function [isBHDR, fstem] = MRIisBHDR(fspec)
% [isBHDR fstem] = MRIisBHDR(fspec)
%
% Determines whether the given file spec is an BHDR
% file based soley on its extension. BHDR is a bshort or
% bfloat.
%
% Returns non-zero if fspec is the name of an BHDR file,
%   Returns 1 if it has a .bhdr extension.
% Returns 0 otherwise.
%


%
% MRIisBHDR.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
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


fstem = [];
isBHDR = [];
if(nargin ~= 1)
  fprintf('[isBHDR fstem] = MRIisBHDR(fspec)\n');
  return;
end

isBHDR = 0;
if(length(fspec) < 6) return; end

ext = fspec(end-4:end);

if(strcmp(ext,'.bhdr')) 
  isBHDR = 1; 
  fstem = fspec(1:end-5);
end

return;













