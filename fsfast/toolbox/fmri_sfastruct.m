function sfa = fmri_sfastruct()
% sfa = fmri_sfastruct()
%


%
% fmri_sfastruct.m
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

sfa.version = 1;
sfa.analysistype = 'average';

sfa.stimtype  = '';
sfa.direction = 'pos';
sfa.ncycles = 0;

sfa.dof = 0;
sfa.TR = 0;
sfa.ntp = 0;
sfa.Trun = 0;
sfa.fundamental = 0;
sfa.delay = 0;
sfa.delay_stem = '';
sfa.freqskip = 3;
sfa.skirtskip = 1;

sfa.isignal  = [];
sfa.inoise   = [];
sfa.iexclude = [];

sfa.nrows   = 0;
sfa.ncols   = 0;
sfa.slice_delay = [];
sfa.firstslice = 0;
sfa.nslices = 0;

sfa.meanval = 0;
sfa.rescale_target = 0;
sfa.rescale_factor = 0;
sfa.nskip = 0;
sfa.hanrad = 0;
sfa.detrend = 0;

sfa.pxform = 'log10';
sfa.stattype = '';      % perharm, sumharm, cumsumharm
sfa.statharmonics = [];

sfa.infiles = [];

return;
