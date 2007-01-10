function slines = readln(fid,nl)
%
% slines = readln(fid,<nlines>)
%
%
%


%
% readln.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
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


F = fread(fid);
lF = length(F);

if( isempty(F) )
  warning('File appears to be empty');
  slines = [];
  return;
end

indNL = find(F == 10); % indicies of new lines
nNL = length(indNL);

if(max(indNL) == lF) % last char is NL
  nLines = nNL;
else
  nLines = nNL + 1;
end

if(nLines == 1)
  slines = setstr(F');
  return;
end

if(nargin == 2) nlmax = nl;
else            nlmax = nLines;
end

if(nlmax > nLines) nlmax = nLines; end

sF = setstr(F');
clear F;
slines = [];
b = 1;
for n = 1:nlmax

  % need this in case file does not end in a newline
  if(n > nNL) e = lF;
  else        e = indNL(n);
  end

  if(e > b + 1)  % skip lines with only NLs
    s = sF(b:e); 
    slines = strvcat(slines,s);
  end

  if(n <= nNL) b = indNL(n) + 1; end

end

return;
