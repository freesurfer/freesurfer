function st = fast_par2st(par,condno,trun)
% st = fast_par2st(par,condno,trun)
%
% Convert a paradigm to  stimulus timing matrix
%
% par - two columns, fully specified timing
%   1. Stimulus onset time
%   2. Condition Number
% 
%
% Col 1: stimulus onset time (sec)
% Col 2: stimulus duration (sec)
% Col 3: stimulus weight (always 1)
%
% $Id: fast_par2st.m,v 1.1 2006/07/19 04:46:29 greve Exp $

st = [];

if(nargin ~= 3)
  fprintf('st = fast_par2st(par,condno,trun)\n');
  return;
end

nstim = size(par,1);

ind = find(par(:,2) == condno);
if(isempty(ind))
  fprintf('ERROR: cannot find condition %d in paradigm\n',condno);
  return;
end
t = par(ind,1);
npres = length(ind);

st = zeros(npres,3);
for nthpres = 1:npres
  tpres = t(nthpres);
  if(ind(nthpres)+1 <= nstim)
    tnext = par(ind(nthpres)+1,1);
  else
    tnext = trun;
  end
  st(nthpres,1) = tpres;
  st(nthpres,2) = tnext-tpres;
  st(nthpres,3) = 1;
end

return;
