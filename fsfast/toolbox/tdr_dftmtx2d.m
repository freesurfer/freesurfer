function [F kA1d kB1d] = tdr_dftmtx2d(arg1,arg2,accel)
% F = tdr_dftmtx2d(mask,<[kAaccel kBaccel]>);
% F = tdr_dftmtx2d(Nrr,Nrc,<[kAaccel kBaccel]>);
% 
% F is Nk-by-Nv matrix
% A is the faster of the two k-space dim, which means
% that if the k-space data should not be permuted when
% read into matlab.
%
% See also ~/links/sg1/propeller2/reconpropss.m
% 


%
% tdr_dftmtx2d.m
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

F = [];
if(nargin < 1 | nargin > 3)
  fprintf('F = tdr_dftmtx2d(mask,   <[kAaccel kBaccel]>)\n');
  fprintf('F = tdr_dftmtx2d(Nrr,Nrc,<[kAaccel kBaccel]>)\n');
  return;
end

if(length(size(arg1)) == 2)
  % Arg1 is a matrix, interpret as a mask
  if(nargin > 2)
    fprintf('ERROR: only on other input arg with mask\n');
    return;
  end
  mask = arg1;
  accel = [1 1];
  if(nargin == 2) accel = arg2; end
  [Nrr Nrc] = size(mask);
else
  % Arg1 is a scalar, intpret as Nrr
  if(nargin < 2)
    fprintf('ERROR: need Nrc with Nrr\n');
    return;
  end
  Nrr = arg1;
  Nrc = arg2;
  if(nargin == 3) accel = arg2; 
  else accel = [1 1];
  end
  mask = ones(Nrr,Nrc);
end
indmask = find(mask);

if(length(accel) ~= 2)
  fprintf('ERROR: accel must be a two component vector\n');
  return;
end
kAaccel = accel(1);
kBaccel = accel(2);

Nkr = Nrr;
Nkc = Nrc; 

% ColNo0 is the column number, centered at Nrc/2+1
ColNo0 = [0:Nrc-1];
ColNo0 = ColNo0 - ColNo0(round(Nrc/2) + 1);

% RowNo0 is the row number, centered at Nrr/2+1
RowNo0 = [0:Nrr-1]';
RowNo0 = RowNo0 - RowNo0(round(Nrr/2) + 1);

% ColNo and RowNo are images of the the col and row numbers at each voxel
ColNo = repmat(ColNo0, [Nrr 1]);  
RowNo = repmat(RowNo0, [1 Nrc]);
ColNo1d = ColNo(indmask)';
RowNo1d = RowNo(indmask)';

% Faster dimension
% kA0 is the fractional row number at each kimage voxel
% Why is this not centered at Nkr/2+1?
kA0 = 2*pi*[0:Nkr-1]/Nkr;
kA0 = kA0 - kA0(round(Nkr/2) + 1);
kA0 = kA0(1:kAaccel:end);
Nkr = length(kA0);
kA0 = reshape(kA0,[Nkr 1]); % column vector

% Slower dimension
% kB0 is the fractional col number at each kimage voxel
% Why is this not centered at Nkc/2+1?
kB0 = 2*pi*[0:Nkc-1]/Nkc;
kB0 = kB0 - kB0(round(Nkc/2) + 1);
kB0 = kB0(1:kBaccel:end);
Nkc = length(kB0);
kB0 = reshape(kB0,[1 Nkc]); % row vector

% kA and kB are images of the the kcol and krow numbers at each voxel
kA = repmat(kA0,[1 Nkc]);
kB = repmat(kB0,[Nkr 1]);
nkv = prod(size(kA)); 
kA1d = reshape(kA,[nkv 1]);
kB1d = reshape(kB,[nkv 1]);

%plot(kB1d,kA1d,'.',kB1d(1),kA1d(1),'r*',kB1d(2),kA1d(2),'g*')
%keyboard

% Compute the phase with an outer product
% kFastest with rFastest
%phi = kA1d*RowNo1d + kB1d*ColNo1d;
% Encoding matrix
%F = exp(-i*phi);

% Do it in one step, may be faster
F = exp(-i*(kA1d*RowNo1d + kB1d*ColNo1d));
 
return;




