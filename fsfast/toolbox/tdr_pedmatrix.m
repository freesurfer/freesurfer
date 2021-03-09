function pedmat = tdr_pedmatrix(TE,echospacing,delsamp,tDwell,nrows,ncols,perev)
% pedmat = tdr_pedmatrix(TE,echospacing,delsamp,tDwell,nrows,ncols,<perev>)
% pedmat = tdr_pedmatrix(epipar);
%
% Computes the EPI post-excitation delay (PED) of each sample in a
% k-space image.  The row/col in the PED matrix corresponds to that of
% the kspace image without any reversals (readout or phase encode)
% applied. Any such operations applied to the kspace image should then
% be applied to the PED matrix. It is assumed that the even lines were
% acquired in reverse readout.
%
% The perev indicates that k-space was traversed in the reverse
% phase encode direction. This changes where the center of k-space is
% in the PED matrix. No reversals or flipping of the PED matrix are
% applied even when the perev is set. The upper left voxel in the
% PED matrix will always have the shortest PED because that was
% the first voxel to be acquired. Setting perev has a subtle
% effect. It only shifts the sample times so that the center of
% k-space is at nrows/2-1 instead of nrows/2+1.
%
% TE - echo time - time at which the center of k-space is traversed.
% echospacing - time between lines/echoes
% delsamp - time to first sample after the start of the ramp
% tDwell - time between ADC samples
% nrows - number of rows in the raw k-space image
% ncols - number of columns in the raw k-space image (includes oversampling)
% perev - k-space was traversed in the reverse direction (PED
%           matrix is not flipped UD). If perev is not given or is
%           empty, no reversal is assumed.
%
%
%


%
% tdr_pedmatrix.m
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

pedmat = [];

if(nargin ~= 1 & nargin ~= 6 & nargin ~= 7)
  fprintf('pedmat = tdr_pedmatrix(TE,echospacing,delsamp,tDwell,nrows,ncols,<perev>)\n');
  fprintf('pedmat = tdr_pedmatrix(epipar)\n');
  return;
end

if(nargin == 1)
  % epipar struct
  epipar = TE; % Just first arg
  TE = epipar.TE;
  echospacing = epipar.echospacing;
  tDwell = epipar.tDwell;
  ncols = epipar.nkcols;
  nrows = epipar.nkrows;
  delsamp = epipar.tDelSamp;
  if(~isfield(epipar,'perev')) epi.perev = 0; end
  perev = epi.perev;
end

if(exist('perev') ~= 1) perev = []; end
if(isempty(perev)) perev = 0; end

% Compute the row where center of k-space is traversed%
if(~perev) % not reversed
  rcenter = nrows/2 + 1;
else % reversed
  rcenter = nrows/2 - 1;
end

% Column center is always here
ccenter = ncols/2 + 1;

% Time to the echo of (ie, center of) the first line
L1CenterPED = TE - (echospacing*(rcenter-1));

% Time to the start of the rampup of the first line
FirstRampPED = L1CenterPED - echospacing/2;

% Time at which each column in a line is sampled wrt start of Ramp
col0PED = tDwell * [0:ncols-1] + delsamp;

% Time at which the ramp of each row starts wrt the start of first Ramp
row0PED = echospacing * [0:nrows-1]';

% Relative time at which each row/col is sampled wrt 
% the start of first Ramp.
pedmat0 = repmat(row0PED,[1 ncols]) + repmat(col0PED,[nrows 1]);

% Reverse readout for even lines
evenrows = 2:2:nrows;
pedmat0(evenrows,:) = fliplr(pedmat0(evenrows,:));

% Relative time at which each row/col is sampled wrt 
% the center of the RF excitation pulse
pedmat = pedmat0 + FirstRampPED;

% Note: pedmat(rcenter,ccenter) should equal TE (+/- tDwell/2)
% pedmat(rcenter,ccenter)-TE

return;
