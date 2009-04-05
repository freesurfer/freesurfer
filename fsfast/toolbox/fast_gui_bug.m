function bug = fast_gui_bug(cspec)
% bug = fast_gui_bug(cspec)
%
% This function determines whether a contrast has the FSFAST GUI
% contrast bug (see
% https://surfer.nmr.mgh.harvard.edu/fswiki/FsFastGuiBug) 
% 
% cspec is just the contrast matfile loaded as a structure, eg,
% cspec = load('analysis/contrast.mat');
%
% Returns:
%  0 - no bug
%  1 - bug but contrast is balanced
%  2 - bug and contrast is unbalanced
%
% $Id: fast_gui_bug.m,v 1.1 2009/04/05 23:22:10 greve Exp $

bug = [];
if(nargin ~= 1)
  fprintf('bug = fast_gui_bug(cspec)\n');
  return;
end

bug = 0; 
if(~isfield(cspec,'CondState'))
  % Lack of CondState means that GUI was NOT used.
  return;
end
if(max(abs(cspec.CondState(:))) == 0)
  % CondState==0 means that GUI was NOT used.
  return;
end  
if(~cspec.CNorm)
  % It was not normalized, so not a problem.
  return;
end  

% Now determine if the WCond is actually correctly normalized by
% normalizing it again. If it is properly normalized, then the
% values should not change.
WCondNorm = fast_norm_con2(cspec.WCond);
m = max(abs(WCondNorm(:)-cspec.WCond(:)));
if(m < .00000001) 
  % Values did not change.
  return; 
end

% If it gets here, then there is a problem.  Determine whether it is a
% balanced or unbalanced design to see how big the problem is.
npos = length(find(cspec.WCond > 0));
nneg = length(find(cspec.WCond < 0));
if( (npos > 0 & nneg > 0) & npos ~= nneg ) 
  % Unbalanced, bad.
  bug = 2;
else 
  % Balanced, not so bad.
  bug = 1;
end


return;

%----------------------------------------------------------------------
function Cnorm = fast_norm_con2(C)
% This is a replication of fast_norm_con.m. It is included here so
% that this file can be self-contained.
%
% Cnorm = fast_norm_con(C)
%
% Make sure that the positives of each row sum to 1 and that the
% negatives of each row sum to -1. This also assures that each row
% sums to zero if there are positives and negatives in the row.
%
% This is a replicated of code found in fast_contrastmtx.m
%


Cnorm = [];
if(nargin ~= 1)
  fprintf('Cnorm = fast_norm_con(C)\n');
  return;
end

Cnorm = C;
for nthrow = 1:size(C,1);
  % positives %
  ind = find(C(nthrow,:)>0);
  if(~isempty(ind))
    xsum = sum(C(nthrow,ind));
    Cnorm(nthrow,ind) = C(nthrow,ind)/xsum;
  end
  % negatives %
  ind = find(C(nthrow,:)<0);
  if(~isempty(ind))
    xsum = sum(C(nthrow,ind));
    Cnorm(nthrow,ind) = C(nthrow,ind)/abs(xsum);
  end
end

return;
