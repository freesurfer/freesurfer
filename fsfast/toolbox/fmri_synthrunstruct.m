function synthrunstruct = fmri_syntrunstruct()

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
