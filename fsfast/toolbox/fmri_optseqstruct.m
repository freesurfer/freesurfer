function oss = fmri_optseqstruct


%
% fmri_optseqstruct.m
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







