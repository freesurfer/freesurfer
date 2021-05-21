function h_is = fmri_havg_is(hs,Misa)
%
% h_is = fmri_havg_is(hs)
% h_is = fmri_havg_is(hs,Misa)
%
% Computes the inter-session HDR average.  hs is
% the HDR for all the sessions over which to average
% (nRows x nCols x Nch x Ns). Misa is the inter-session 
% averaging matricies (Nch x Nch x Ns). If Misa is not
% given, then the inter-session average at each voxel
% at each dealy-point is just the average across all
% the sessions.  h_is will have dimension (nRows x NCols x
% Nch).
%
% See also fmri_avgmtx_is.
%
%


%
% fmri_havg_is.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'h_is = fmri_havg_is(hs,<Misa>)';
  qoe(msg);error(msg);
end


[nRows nCols Nch Ns] = size(hs);

% Misa is not given, just compute the average %
if(nargin == 1)
  h_is = mean(hs,4);
  return;
end

% Only reaches this point if Misa is given %

if(size(Misa,1) ~= Nch | size(Misa,3) ~= Ns)
  fprintf('size(Misa) = '); fprintf(1,'%d  ',size(Misa));fprintf(1,'\n');
  msg = sprintf('Misa size is inconsistent with Nch=%d or Ns=%d',Nch,Ns);
  qoe(msg);error(msg);
end

nVoxels = nRows*nCols;
hs = reshape(hs,[nVoxels Nch Ns]);
hs = permute(hs,[2 1 3]);

h_is = zeros(Nch,nVoxels);
for s = 1:Ns,
  h_is = h_is + Misa(:,:,s)*hs(:,:,s);
end

h_is = reshape(h_is',[nRows nCols Nch]);

return;
