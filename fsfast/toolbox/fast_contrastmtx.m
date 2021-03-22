function R = fast_contrastmtx(TER,TW,TPS,nConds,SumConds,WConds,SumDelays,WDelays,RmPrestim,CNorm)
% 
% R = fast_contrastmtx(TER,TW,TPS,nConds,SumConds,WConds,
%                      SumDelays,WDelays,RmPrestim,CNorm)
%
% R = fast_contrastmtx(cspec);
%
% SumConds = 1, forces Conditions to be weighted by WConds and summed.
%
% WConds - condition weighting vector. Ignored if SumConds = 0. If
%  [], replaced with ones (ie, forces a simple average).
%
% SumDelays = 1, forces Delays to be weighted by WDelays and summed.
%   Forces the contrast matrix to be a vector.
%
% WDelays - delay weighting vector. Ignored if SumDelays = 0. If
%  [], replaced with ones (ie, forces a simple average).
%
% RmPrestim = 1 subtract prestimulus average; will replace
%   the prestim components of WDelays.
%
% CNorm = 1: set weights such that each row sums to 0
%         0: do not change from the original


%
% fast_contrastmtx.m
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

R = [];

if(nargin ~= 10 & nargin ~= 1)
  fprintf('R = fast_contrastmtx(TER,TW,TPS,nConds,SumConds,WConds,SumDelays,WDelays,RmPrestim,CNorm)');
  fprintf('R = fast_contrastmtx(cspec);\n');
  return;
end

if(nargin == 1)
  cspec     = TER;
  TER       = cspec.TER;
  TW        = cspec.TimeWindow;
  TPS       = cspec.TPreStim;
  nConds    = cspec.NCond;
  SumConds  = cspec.sumconds;
  WConds    = cspec.WCond;
  SumDelays = cspec.sumdelays;
  WDelays   = cspec.WDelay;
  RmPrestim = cspec.RmPreStim;
  CNorm     = cspec.CNorm;
end


nDelays = round(TW/TER);

if(SumConds)
  if(isempty(WConds)) WConds = ones(nConds,1); end
  if(length(WConds) ~= nConds) 
     fprintf('ERROR: WConds length = %d, should be %d\n',...
	     length(WConds),nConds);
     return;
  end
end

R = [];
for nthCond = 1:nConds

  RCond = fast_condctrstmtx(TER,TW,TPS,SumDelays,WDelays,RmPrestim);

  if(isempty(RCond)) R = []; return; end

  if(SumConds)
    R = [R WConds(nthCond)*RCond];
  else
    if(WConds(nthCond) ~= 0)
      R0       = zeros(size(RCond));
      Rprepad  = repmat(R0,[1 (nthCond-1)]);
      Rpostpad = repmat(R0,[1 (nConds-nthCond)]);
      Rc = [Rprepad WConds(nthCond)*RCond Rpostpad];
      R = [R; Rc];
    end
  end
end


if(CNorm==1)
  % Make sure that the positives of each row sum to 1
  % and that the negatives of each row sum to -1. This
  % also assures that each row sums to zero if there
  % are positives and negatives in the row.
  fprintf('Normalizing contrast matrix\n');
  % This code is replicated in fast_norm_con.m
  for nthrow = 1:size(R,1);
    % positives %
    ind = find(R(nthrow,:)>0);
    if(~isempty(ind))
      xsum = sum(R(nthrow,ind));
      R(nthrow,ind) = R(nthrow,ind)/xsum;
    end
    % negatives %
    ind = find(R(nthrow,:)<0);
    if(~isempty(ind))
      xsum = sum(R(nthrow,ind));
      R(nthrow,ind) = R(nthrow,ind)/abs(xsum);
    end
  end
end

return;

