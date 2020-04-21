function irepistruct = irepitiming(irepistruct)
% irepistruct = irepitiming(irepistruct)
% 
% See irepistructure.m for more info on the struct parameters
%

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
  if(s.Slice1PreInv==0)
    % Inversion
    t = t - s.PreInv;    
    t = t + s.InvDur; % duration of inversion pulse
    tInv = t;
    s.tEvent = [s.tEvent; tInv];
    s.FlipEvent = [s.FlipEvent; InvFlip];
    s.IsReadOut = [s.IsReadOut; 0];
    s.TI = [s.TI; 0];
    s.AcqSliceNo = [s.AcqSliceNo; 0];
    s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
    t = t + s.TI1;
  end
  % Readout (but no actual acq because these are dummies)
  nthslice = 0;
  for sliceno = SliceOrderD
    nthslice = nthslice + 1;
    s.tEvent = [s.tEvent; t];
    s.EventSliceNo = [s.EventSliceNo; sliceno];
    s.FlipEvent = [s.FlipEvent; s.ROFlip];
    s.IsReadOut = [s.IsReadOut; 0]; % 0 = no acq
    s.TI = [s.TI; 0];
    s.AcqSliceNo = [s.AcqSliceNo; 0];
    t = t + s.TBS;
    if(nthslice == 1 & s.Slice1PreInv)
      % Inversion
      t = t - s.PreInv;
      t = t + s.InvDur; % duration of inversion pulse
      tInv = t;
      s.tEvent = [s.tEvent; tInv];
      s.FlipEvent = [s.FlipEvent; InvFlip];
      s.IsReadOut = [s.IsReadOut; 0];
      s.TI = [s.TI; 0];
      s.AcqSliceNo = [s.AcqSliceNo; 0];
      s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
      t = t + s.TI1;
    end
  end
end

% Non-dummies
for nthtp = 1:s.ntp
  if(s.Slice1PreInv==0)
    % Inversion before all slices
    t = t - s.PreInv;
    t = t + s.InvDur; % duration of inversion pulse
    tInv = t;
    s.tEvent = [s.tEvent; tInv];
    s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
    s.FlipEvent = [s.FlipEvent; InvFlip];
    s.IsReadOut = [s.IsReadOut; 0];
    s.TI = [s.TI; 0];
    s.AcqSliceNo = [s.AcqSliceNo; 0];
    t = t + s.TI1;
  end
  % Read out
  SliceOrderR = circshift([1:s.nslices]',-(nthtp-1)*s.skip)';
  nthslice = 0;
  for sliceno = SliceOrderR
    nthslice = nthslice + 1;
    s.tEvent = [s.tEvent; t];
    s.EventSliceNo = [s.EventSliceNo; sliceno];
    s.FlipEvent = [s.FlipEvent; s.ROFlip];
    s.IsReadOut = [s.IsReadOut; 1];
    s.TI = [s.TI; t-tInv];
    s.AcqSliceNo = [s.AcqSliceNo; nthslice];
    t = t + s.TBS;
    if(nthslice == 1 & s.Slice1PreInv)
      % Inversion, when 1st slice occurs before inversion
      t = t - s.PreInv;
      t = t + s.InvDur; 
      tInv = t;
      s.tEvent = [s.tEvent; tInv];
      s.EventSliceNo = [s.EventSliceNo; -1]; % Applies to all slices
      s.FlipEvent = [s.FlipEvent; InvFlip];
      s.IsReadOut = [s.IsReadOut; 0];
      s.TI = [s.TI; 0];
      s.AcqSliceNo = [s.AcqSliceNo; 0];
      t = t + s.TI1;
    end
  end % slice
end

irepistruct = s;

return;

