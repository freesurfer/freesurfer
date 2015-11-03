function irepistruct = irepitiming(irepistruct)
% irepistruct = irepitiming(irepistruct)
% 
% See irepistructure.m for more info on the struct parameters
%
% $Id: irepitiming.m,v 1.1 2015/10/30 21:59:56 greve Exp $

InvFlip = 180;

s = irepistruct; % copy into s for easy handling
s.tEvent = [];
s.EventSliceNo = []; % spatial/anatomical slice no
s.AcqSliceNo = []; % after inversion
s.FlipEvent = [];
s.IsReadOut = [];

t = 0;

% Dummies.
SliceOrderD = [1:s.nslices]; % Dummies do not permute
for nthtp = 1:s.ndummies
  t = t + s.InvDur; % Prefill = duration of inversion pulse
  s.tEvent = [s.tEvent; t];
  % Inversion
  s.FlipEvent = [s.FlipEvent; InvFlip];
  s.IsReadOut = [s.IsReadOut; 0];
  s.AcqSliceNo = [s.AcqSliceNo; 0];
  s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
  t = t + s.TI1;
  % Readout (but no actual acq because these are dummies)
  nthslice = 0;
  for sliceno = SliceOrderD
    nthslice = nthslice + 1;
    s.tEvent = [s.tEvent; t];
    s.EventSliceNo = [s.EventSliceNo; sliceno];
    s.FlipEvent = [s.FlipEvent; s.ROFlip];
    s.IsReadOut = [s.IsReadOut; 0]; % 0 = no acq
    s.AcqSliceNo = [s.AcqSliceNo; 0];
    t = t + s.TBS;
  end
end

for nthtp = 1:s.ntp
  t = t + s.InvDur; % Prefill = duration of inversion pulse
  s.tEvent = [s.tEvent; t];
  % Inversion
  s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
  s.FlipEvent = [s.FlipEvent; InvFlip];
  s.IsReadOut = [s.IsReadOut; 0];
  s.AcqSliceNo = [s.AcqSliceNo; 0];
  t = t + s.TI1;
  % Read out
  SliceOrderR = circshift([1:s.nslices]',-(nthtp-1)*s.skip)';
  nthslice = 0;
  for sliceno = SliceOrderR
    nthslice = nthslice + 1;
    s.tEvent = [s.tEvent; t];
    s.EventSliceNo = [s.EventSliceNo; sliceno];
    s.FlipEvent = [s.FlipEvent; s.ROFlip];
    s.IsReadOut = [s.IsReadOut; 1];
    s.AcqSliceNo = [s.AcqSliceNo; nthslice];
    t = t + s.TBS;
  end
end

irepistruct = s;

return;

