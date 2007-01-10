function st = fast_par2st(par,condno,trun,stimdur)
% st = fast_par2st(par,condno,trun,<stimdur>)
%
% Convert a paradigm to stimulus timing matrix
%
% par - two columns
%   1. Stimulus onset time
%   2. Condition Number
% 
% st is a matrix with:
%  Col 1: stimulus onset time (sec)
%  Col 2: stimulus duration (sec)
%  Col 3: stimulus weight (always 1)
% If condno is not found in par, then st=[]
%
% If stimdur is specified, then all stimuli will have that as the
% stimulus duration (column 3)
%
%


%
% fast_par2st.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.3 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

st = [];

if(nargin ~= 3 & nargin ~= 4)
  fprintf('st = fast_par2st(par,condno,trun,<stimdur>)\n');
  return;
end

nstim = size(par,1);

ind = find(par(:,2) == condno);
if(isempty(ind))
  %fprintf('ERROR: cannot find condition %d in paradigm\n',condno);
  % Don't print an error here, in case allowing some conditions to 
  % not exist in some runs. Just return an empty.
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

% If stimdur specified, use with all presentations
if(nargin == 4) st(:,2) = stimdur; end

return;
