function sfa = fmri_sfastruct()
% sfa = fmri_sfastruct()
% $Id: fmri_sfastruct.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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
