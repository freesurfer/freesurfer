function Q = vqm(abc,voxsize,voldim,base)
% 
% Q = vqm(abc,voxsize,voldim,<base>)
%
% Volume quantization matrix. Q computes the x,y,z of a voxel given 
% the subscript of the  voxel in the volume.
%
% abc - 3 character string indicating which physical dimension maps to
%       the  fast (a), medium (b), and slow (c) indicies. The string also 
%       indicates the sign of the derivative of the physical dimension 
%       with respect to its index. Valid characters are (case insensitive):
%          R - x dimension, dx/di > 0 (index increases Left to Right)
%          L - x dimension, dx/di < 0 (index increases Right to Left)
%          A - y dimension, dy/di > 0 (index increases Post to Ant)
%          P - y dimension, dy/di > 0 (index increases Ant to Post)
%          S - z dimension, dz/di > 0 (index increases Infer to Super)
%          I - z dimension, dz/di > 0 (index increases Super to Infer)
%
% voxsize - 3 element vector of the (absolute) physical length of
%       a voxel side.  Each element corresponds to the phyical 
%       dimension in the abc string.
%
% voxdim - 3 element vector of the number of voxels in each dimension.
%       Each element corresponds to the phyical dimension in the abc string.
%
% base - the value of the first index (default is 0).
%
% Example, default for MGH-NMR FreeSsurfer Corronals:
%
%     Q = vqm('LIA', [1 1 1], [256 256 256], 0);
% 
% Example, default for MGH-NMR FreeSsurfer Functionals
%     [R subj inres betres] = fmri_readreg('register.dat'); 
%     Qf = vqm('LIA', [inres inres betres], [64 64 24], 0);
%
%


%
% vqm.m
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

if(nargin ~= 3 & nargin ~= 4)
  msg = 'USAGE: Q = vqm(abc,voxsize,voldim,<base>)';
  error(msg);
end

if( length(voxsize) ~= 3)
  msg = sprintf('voxsize dimension is wrong (%d ~= 3)',length(voxsize));
  error(msg);
end

[ok msg] = checkabc(abc);
if(~ok)  error(msg); end

if(nargin == 3) base = 0; end

voxsize = abs(reshape1d(voxsize));
voxdim  = reshape1d(voldim);
voxsign = abc2voxsign(abc);

Q = zeros(4,4);
Q(4,4) = 1;

for n = 1:3,
  s = abc(n);
  switch upper(s)
    case {'R','L'}, v = [1 0 0]; m = 1;
    case {'P','A'}, v = [0 1 0]; m = 2;
    case {'I','S'}, v = [0 0 1]; m = 3;
  end

  f = voxsign(m) * voxsize(m) ;
  Q(n,[1:3]) = f*v;
  Q(n,4) = -voxsign(m) * (voxsize(m)*(voxdim(m)+base)/2);
end

return;
%---------------------------------------%
%------------>>>>>\./<<<<<<-------------%
%---------------------------------------%


%---------------------------------------%
function voxsign = abc2voxsign(abc)

  voxsign = zeros(3,1);

  for n = 1:3,
    s = abc(n);

    switch upper(s)
      case {'R'}, voxsign(n) = +1;
      case {'A'}, voxsign(n) = +1;
      case {'S'}, voxsign(n) = +1;
      case {'L'}, voxsign(n) = -1;
      case {'P'}, voxsign(n) = -1;
      case {'I'}, voxsign(n) = -1;
      otherwise,
        msg = sprintf('Unknown code in %s',abc);
        error(msg);
    end

  end
return;
%---------------------------------------%

%---------------------------------------%
function [ok,msg] = checkabc(abc)

  ok = 1;
  msg = [];
  abc0 = abc;
  abc = upper(abc);
  for n = 1:3,

    v = [1:3];
    i = find(v ~= n);
    v = v(i);

    switch abc(n)
      case {'R','L'}, 
        i = find(abc(v)=='R');
        j = find(abc(v)=='L');
      case {'A','P'}, 
        i = find(abc(v)=='A');
        j = find(abc(v)=='P');
      case {'S','I'}, 
        i = find(abc(v)=='S');
        j = find(abc(v)=='I');
      otherwise,
        ok = 0;
        msg = sprintf('Unknown code abc in %s',abc0);
        return;
    end

   if( ~isempty(i) | ~isempty(j) )
     ok = 0;
     msg = sprintf('Duplicate abc code in %s',abc0);
     return;
   end

  end

return


%---------------------------------------%
%----- Reshape into a column vector ----%
%---------------------------------------%
function y = reshape1d(x)
  y = reshape(x, [prod(size(x)) 1]);
return;
%---------------------------------------%
