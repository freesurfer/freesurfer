function hdrdat = fmri_hdrdatstruct()
% hdrdat = fmri_hdrdatstruct()
%


%
% fmri_hdrdatstruct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.3 $
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

hdrdat.Version = 3;
hdrdat.TR = 0;
hdrdat.TER = 0;
hdrdat.TimeWindow = 0;
hdrdat.TPreStim = 0;
hdrdat.Nc = 0; % Total number of conditions, incl fix
hdrdat.Nh = 0;
hdrdat.Nnnc = 0;
hdrdat.DOF  = 0;
hdrdat.Npercond = [];
hdrdat.Nruns = 0;
hdrdat.Ntp = 0;
hdrdat.Nrows = 0;
hdrdat.Ncols = 0;
hdrdat.Nskip = 0;
hdrdat.DTOrder = 0;
hdrdat.RescaleFactor = 1;
hdrdat.HanningRadiues = 0;
hdrdat.BrainAirSeg = 1;
hdrdat.GammaFit = 0;
hdrdat.gfDelta = [];
hdrdat.gfTau = [];
hdrdat.NullCondId = 0;
hdrdat.SumXtX = [];
hdrdat.nNoiseAC = 0;
hdrdat.hCovMtx = [];
hdrdat.CondIdMap = [];

hdrdat.LPFFlag = 0;   % 1 = use global low-pass filter
hdrdat.HPF = [-1 -1]; % global high-pass using inverse of ar1w model
hdrdat.WhitenFlag = 0;    % per-run high-pass

% FS-FAST-related parameters 
hdrdat.runlist = [];
hdrdat.funcstem = [];
hdrdat.parname = [];
hdrdat.extregstem = [];
hdrdat.nextreg = [];
hdrdat.extregortho = [];


return;

