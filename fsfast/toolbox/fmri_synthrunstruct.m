function synthrunstruct = fmri_syntrunstruct()


%
% fmri_synthrunstruct.m
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

srs.TR         = 2;
srs.Rss        = 1;
srs.TER        = srs.TR;
srs.TimeWindow = 20;
srs.delta      = [2.25];
srs.tau        = [1.25];
srs.SNR        = -1.0; % infinite
srs.Offset     = 1000;
srs.Trend      = 0;
srs.PercentSigChange = 5; % Percent of Offset %
srs.RnnMaxDelay  = 20.0; % seconds
srs.alpha = .75;
srs.rho   = .88;
srs.Nrows = 64;
srs.Ncols = 64;
srs.Ntp   = 128;
srs.Nskip = 0;
srs.Nts   = srs.Ntp*srs.Rss;
srs.Nslices    = 1;
srs.ROI        = [];
srs.NPerCond   = [64 64];
srs.traceiXtX  = 0;
srs.Par        = [];
srs.Nsearch    = 100;
srs.Nswap      = -1;
srs.Sig = [];
srs.SigMean = 0;
srs.SigVar  = 0;
srs.SigSS  = 0;
srs.Seed   = -1;

synthrunstruct = srs;

return;
