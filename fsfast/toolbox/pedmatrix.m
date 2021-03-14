function pedmat = pedmatrix(echo1ped,echospacing,delsamp,tDwell,nrows,ncols,nechoes)
% pedmat = pedmatrix(echo1ped,echospacing,delsamp,tDwell,nrows,ncols,nechoes) 
%
% Computes the post-excitation delay of each sample in a k-space
% image, including multiple echos. The k-space image is
% nrows-by-ncols with even lines going in the same direction as odd
% lines (ie, if the even lines were reversed during acquisition,
% then they have been reversed again in the k-space image). If
% there are multiple echos, then it is assumed that the even echos
% were collected in reverse order (but have been reversed again in
% the kspace image).
%
% echo1ped - time to the first echo of the first line
% echospacing - time between lines/echoes
% delsamp - time to first sample after the start of the ramp
% tDwell - time between ADCs
% 
%
%


%
% pedmatrix.m
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

if(nargin ~= 7)
  fprintf('pedmat = pedmatrix(echo1ped,echospacing,delsamp,tDwell,nrows,ncols,nechoes)\n');
  return;
end

% Time from the RF Pulse to the start of very first rampup
FirstRampPED = echo1ped - echospacing/2;

% Time at which each column in a line is sampled wrt start of Ramp
tcol0 = tDwell * [0:ncols-1] + delsamp;

% Time at which the ramp of each row starts wrt the start of first Ramp
trow0 = echospacing * [0:nrows-1]';

% Time to traverse k-space for one image
Timage = echospacing * nrows;

% Time of first ramp of each kimage
techo0 = [0:nechoes-1]*Timage;

% Time of first and second major echoes
% techo1 = trow0(nrows/2+1) + echo1ped;
%   Should be at pedmat(nrows/2+1,ncols/2+1,1)
% techo2 = techo1 + echospacing*(nrows-1); % nrows-1 for rev readout
%   Should be at pedmat(nrows/2+1,ncols/2+1,2)

% Relative time at which each row/col is sampled
pedmat0 = repmat(trow0,[1 ncols]) + repmat(tcol0,[nrows 1]);

% Reverse Readout
evenrows = [2:2:nrows];
pedmat0(evenrows,:) = fliplr(pedmat0(evenrows,:));

for echo = 1:nechoes
  if(rem(echo,2) ~= 0) 
    % odd numbered echo
    pedmat(:,:,echo) = pedmat0;
  else
    % even numbered echo - reverse row order
    pedmat(:,:,echo) = flipud(pedmat0);
  end
  pedmat(:,:,echo) = pedmat(:,:,echo) + techo0(echo);
end

pedmat = pedmat + FirstRampPED;

return;
