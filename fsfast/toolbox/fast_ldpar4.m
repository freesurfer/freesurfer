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
% a. Two cols:   onset time, condition number
% b. Three cols: onset time, condition number, duration
% c. Four cols:  onset time, condition number, duration, stringid
% d. Four cols:  onset time, condition number, duration, weight
% e. Three cols: onset time, condition number, duration
% f. Three cols: onset time, condition number, stringid

%
% fast_ldpar4
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/06/24 19:07:04 $
%    $Revision: 1.5 $
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

nthrow = 0;
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
  if(count < 2)
    fprintf('ERROR: %s is not correctly formatted. ',par4file);
    fprintf('Line %d only has %d items\n',nthrow+1,count);
    fclose(fp);
    par4 = [];
    return;
  end
  if(nthrow == 0) prevcount = count; end
  if(count ~= prevcount & (count < 5 & prevcount < 5))
    % It is ok for different lines have different numbers of cols
    % as long as they have atleast 4 columns. This should prevent
    % the situation where you have 3 cols in one row and 4 in another
    fprintf('ERROR: %s is not correctly formatted. ',par4file);
    fprintf('Line %d has %d items, 1st line had %d\n',...
	    nthrow+1,count,prevcount);
    fclose(fp);
    par4 = [];
    return;
  end
  nthrow = nthrow + 1;
  
  % First column - onset time
  item = sscanfitem(tline,1);
  tonset = sscanf(item,'%f');
  par4(nthrow,1) = tonset;
  
  % Second column - stimulus id
  item = sscanfitem(tline,2);
  condid = sscanf(item,'%d');
  par4(nthrow,2) = condid;

  if(count < 3) continue; end % Only two cols
  
  % Third column
  item = sscanfitem(tline,3);
  duration = sscanf(item,'%f');
  if(isempty(duration)) continue; end % 3rd col is a string
  par4(nthrow,3) = duration;
  
  if(count < 4) continue; end % Only three cols
  
  % Fourth column
  item = sscanfitem(tline,4);
  weight = sscanf(item,'%f');
  if(isempty(weight)) continue; end % 4th col is a string
  % If it gets here, assume it is the weight
  par4(nthrow,4) = weight;

end % while (1)

fclose(fp);

% optseq used to make the last duration a large negative number
if(size(par4,2) > 2 )
  duration = par4(:,3);
  indneg = find(duration < 0);
  if(length(indneg) > 0)
    if(length(indneg) > 1 | indneg(1) ~= size(par4,1))
      fprintf('ERROR: parfile %s has multiple rows with negative duration.\n',...
	      par4file);
      par4 = [];
      return;
    end
    fprintf('WARNING: the duration in the last row in\n');
    fprintf('parfile %s has a negative duration.\n',par4file);
    fprintf('I am going to interpret this as a two col par\n');
    par4 = par4(:,1:2);
  end
end

if(size(par4,2) == 2 )
  fprintf('Loaded %s as a par2 file.\n',par4file);
elseif(size(par4,2) == 3 )
  fprintf('Loaded %s as a par3 file.\n',par4file);
else  
  fprintf('Loaded %s as a par4 file.\n',par4file);
end

return;


