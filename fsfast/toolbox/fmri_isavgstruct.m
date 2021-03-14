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
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


