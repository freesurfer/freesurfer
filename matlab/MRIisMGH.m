function isMGH = MRIisMGH(fspec)
% isMGH = MRIisMGH(fspec)
%
% Determines whether the given file spec is an MGH/MGZ
% file based soley on its extension
%
% Returns non-zero if fspec is the name of an MGH file,
%   Returns 1 if it has a .mgh extension.
%   Returns 2 if it has a .mgz extension.
% Returns 0 otherwise.
%


%
% MRIisMGH.m
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


isMGH = [];
if(nargin ~= 1)
  fprintf('isMGH = MRIisMGH(fspec)\n');
  return;
end

isMGH = 0;
if(length(fspec) < 5) return; end

ext = fspec(end-3:end);

if(strcmp(ext,'.mgh')) isMGH = 1; end
if(strcmp(ext,'.mgz')) isMGH = 2; end

return;













