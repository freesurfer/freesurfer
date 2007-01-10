function oss = fmri_optseqstruct


%
% fmri_optseqstruct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

oss = struct(...
  'Nsessions',    1, ... 
  'TR',           0, ... 
  'Ntp',          0, ... % Number of scans per run
  'Trun',         0, ... % Run Duration (seconds)
  'Nruns',        0, ...
  'TER',          0, ...
  'TimeWindow',   0, ...
  'TPreStim',     0, ...
  'TPreScan',     0, ...
  'Nc',           0, ...  % Total number of conditions, incl fix
  'Nnnc',         0, ...  % Number of non-null conditions
  'Nh',           0, ...
  'Npercond',     [], ...
  'Tpercond',     [], ...
  'Tres',         0, ...  % Temporal resolution of the stim onset.
  'Nsearch',      0, ...
  'NMinsearch',   0, ...  % Search at least this many regardless
  'Tsearch',      0, ...
  'PctSearch',    0, ...  % Search until exceeding 
  'DOF',          0,...
  'pforder',      0,... % polynomial trend order
  'GammaParams', [],... % delta and tau
  'FindWorst',    0,... % for testing purposes
  'MaxEffLimit',  inf); % for testing purposes

return;







