function [T2s, fidhat] = tdr_fidt2star(fid,tfid,nfit)
%
% [T2s, fidhat] = tdr_fidt2star(fid,tfid,<nfit>)
%
% Compute a T2star map from the Free Induction Decay (fid)
% map by fitting the first nfit components of log(fid) 
% to a linear model.
%
% $Id: tdr_fidt2star.m,v 1.1 2003/10/20 22:14:43 greve Exp $
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