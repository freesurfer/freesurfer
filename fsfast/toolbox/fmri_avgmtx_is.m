function Misa = fmri_avgmtx_is(Nch, Ns, XtXs)
%
% Misa = fmri_avgmtx_is(Nch, Ns)
% Misa = fmri_avgmtx_is(Nch, Ns, XtXs)
%
% Compute inter-session averaging matrix.
% Nch - number conditions * number of estimates per event
% Ns  - number of sessions
% XtXs  - list of (sum(x'x)) (one for each session).
%
% Misa will have dimension Nch x Nch x Ns
%
% $Id: fmri_avgmtx_is.m,v 1.1 2003/03/04 20:47:39 greve Exp $

if(nargin ~= 2 & nargin ~= 3)
  msg = 'Misa = fmri_avgmtx_is(Nch, Ns, <XtXs>)';
  qoe(msg);error(msg);
end

Misa = repmat(eye(Nch), [1 1 Ns]);


if(nargin == 2)
  Misa = Misa/Ns;
  return;
end

% Only reaches this point if nargin == 3 %

if(size(XtXs,1) ~= Nch | size(XtXs,3) ~= Ns)
  fprintf('size(XtXs) = '); fprintf(1,'%d  ',size(XtXs));fprintf(1,'\n');
  msg = sprintf('XtXs size is inconsistent with Nch=%d or Ns=%d',Nch,Ns);
  qoe(msg);error(msg);
end

dXtXssum = 0;
for s = 1:Ns,
  dXtXsum = dXtXsum + diag(XtXs(:,:,s));
end

for s = 1:Ns,
  d = diag(XtXs(:,:,s))./dXtXsum;
  Misa(:,:,s) = diag(d);
end

return;
