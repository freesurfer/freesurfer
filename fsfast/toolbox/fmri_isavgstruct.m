function isas = fmri_isavgstruct()
% isas = fmri_isavgstruct()
%
% Creates a structure for inter-subject averaging.
%
%
%


%
% fmri_isavgstruct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:06 $
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


isas = struct(...
  'hdrstem', [], ...
  'hdrgroupid', 0, ...
  'sigfile', [], ...
  'sigformat', 'lnp', ...
  'havg',  [], ...
  'p', [],...
  'w', [] );


