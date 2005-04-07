function [F, rx1d, ry1d, kx1d, ky1d] = tdr_dftmtx2d(Nrr,Nrc,Nkr,Nkc)
% F = tdr_dftmtx2d(Nrr,Nrc,<Nkr>,<Nkc>)
% 
% F is Nv-by-Nv matrix
% See also ~/links/sg1/propeller2/reconpropss.m
% 

if(nargin < 2 | nargin > 4)
  fprrintf('F = tdr_dftmtx2d(Nrr,Nrc,<Nkr>,<Nkc>)\n');
  return;
end

if(~exist('Nkr','var')) Nkr = []; end
if(isempty(Nkr)) Nkr = Nrr; end

if(~exist('Nkc','var')) Nkc = []; end
if(isempty(Nkc)) Nkc = Nrc; end

% rx0 is the column number, centered at Nrc/2+1
rx0 = [0:Nrc-1];
rx0 = rx0 - rx0(round(Nrc/2) + 1);

% ry0 is the row number, centered at Nrc/2+1
ry0 = [0:Nrr-1]';
ry0 = ry0 - ry0(round(Nrr/2) + 1);

% rx and ry are images of the the col and row numbers at each voxel
rx = repmat(rx0, [Nrr 1]);  
ry = repmat(ry0, [1 Nrc]);
nrv = Nrc*Nrr;
rx1d = reshape(rx,[1 nrv]);
ry1d = reshape(ry,[1 nrv]);

% kx0 is the fractional column number at each kimage voxel
% Why is this not centered at Nkc/2+1?
kx0 = 2*pi*[0:Nkc-1]/Nkc;
kx0 = kx0 - kx0(round(Nkc/2) + 1);

% ky0 is the fractional row number at each kimage voxel
% Why is this not centered at Nkr/2+1?
ky0 = 2*pi*[0:Nkr-1]/Nkr;
ky0 = ky0 - ky0(round(Nkr/2) + 1);

% kx and ky are images of the the kcol and krow numbers at each voxel
kx = repmat(kx0,[Nkr 1]);
ky = repmat(ky0,[1 Nkc]);
nkv = Nkc*Nkr;
kx1d   = reshape(kx,[nkv 1]);
ky1d   = reshape(ky,[nkv 1]);
 
% Compute the phase with an outer product
phi = kx1d*rx1d + ky1d*ry1d;

% Encoding matrix
F = exp(-i*phi);
 
return;




