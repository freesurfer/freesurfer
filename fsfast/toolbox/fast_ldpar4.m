function par4 = fast_ldpar4(par4file)
% par4 = fast_ldpar4(par4file)
%
% Load par4-formatted paradigm file
%
% Par4File format is as follows:
%   Column 1: stimulus onset time 
%   Column 2: stimulus condition number
%   Column 3: stimulus duration
%   Column 4: stimulus weight
%   Coumns 5+ are ignored
% Blank lines are ok
% Lines that begin with # are ignormed (comments)
%  
% par4 dimensionality: 4 x nPresentations
%
%

%
% fast_ldpar4
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/06/07 22:12:09 $
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

par4 = [];

if(nargin == 0)
  msg = 'par4 = fast_ldpar4(par4file)';
  qoe(msg);
  error(msg);
end

[fp msg] = fopen(par4file,'r');
if(fp == -1)
  fprintf('ERROR: opening %s\n',par4file);
  fprintf('%s\n',msg);
  return;
end

nthrow = 1;
while(1)
  % scroll through any blank lines or comments 
  % comments are lines that begin with #
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Get count of number of items in the line
  [items count] = sscanf(tline,'%s');
  if(count < 4)
    fprintf('ERROR: %s is not correctly formatted. ',par4file);
    fprintf('Line %d only has %d items\n',nthrow,count);
    fclose(fp);
    par4 = [];
    return;
  end

  for c = 1:4
    item = sscanfitem(tline,c);
    tmp = sscanf(item,'%f');
    if(isempty(tmp))
      fprintf('ERROR: %s is not correctly formatted. ',par4file);
      fprintf('Line %d, item %d appears to be a string\n',nthrow,c);
      fprintf('%s\n',tline);
      fclose(fp);
      par4 = [];
      return;
    end
    par4(nthrow,c) = tmp;
  end

  nthrow = nthrow + 1;
  
end % while (1)

fclose(fp);



return;


