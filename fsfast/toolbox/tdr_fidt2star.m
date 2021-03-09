function [T2s, fidhat] = tdr_fidt2star(fid,tfid,nfit)
%
% [T2s, fidhat] = tdr_fidt2star(fid,tfid,<nfit>)
%
% Compute a T2star map from the Free Induction Decay (fid)
% map by fitting the first nfit components of log(fid) 
% to a linear model.
%
%
%


%
% tdr_fidt2star.m
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

T2s = [];
fidhat = [];

if(nargin ~= 2 & nargin ~= 3)
  fprintf('[T2s, fidhat] = tdr_fidt2star(fid,tfid,<nfit>)\n');
  return;
end

szfid = size(fid);
nechoes = szfid(end);
nv = prod(szfid(1:end-1));

if(exist('nfit') ~= 1) nfit = []; end
if(isempty(nfit)) nfit = nechoes; end
if(nfit > nechoes)
  fprintf('ERROR: nfit = %d > nechoes = %d\n',nfit,nechoes);
  return;
end
if(length(tfid) ~= nechoes)
  fprintf('ERROR: length(tfit) = %d != nechoes = %d\n',...
	  length(tfid),nechoes);
  return;
end

fid = reshape(fid,[nv nechoes])';
X = [ones(nfit,1) -tfid(1:nfit)];
beta = (inv(X'*X)*X')*log(fid(1:nfit,:));
T2s = abs(1./beta(2,:));
T2s = reshape(T2s,[szfid(1:end-1)]);

if(nargout == 2)
  X = [ones(nechoes,1) -tfid(1:nechoes)];
  fidhat = reshape(exp(X*beta)',szfid);
end

return;
