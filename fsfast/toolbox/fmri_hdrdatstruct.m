function hdrdat = fmri_hdrdatstruct()
% hdrdat = fmri_hdrdatstruct()
%


%
% fmri_hdrdatstruct.m
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
hdrdat.HanningRadius = 0;
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

