function [l] = read_label(sname, lname)
% l = read_label(<sname>, lname)
%
% reads the label file 'lname' from the subject 'sname' 
% in the subject's label directory into the vector l
% l will be nvertices-by-5, where each column means:
% (1) vertex number, (2-4) xyz at each vertex, (5) stat
%
% IMPORTANT: the vertex number is 0-based.
% 


%
% read_label.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: fischl $
%    $Date: 2008/02/18 14:37:57 $
%    $Revision: 1.6 $
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


l = [];

if(nargin ~= 2)
  fprintf('l = read_label(<sname>, lname)\n');
  return;
end

if(~isempty(sname))
  sdir = getenv('SUBJECTS_DIR') ;
  fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname) ;
else
	ind = findstr(lname, '.label') ;
	if (length(ind) > 0)
        fname = lname ;
    else
		fname = sprintf('%s.label', lname);
	end
end

% open it as an ascii file
fid = fopen(fname, 'r') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fgets(fid) ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;

fclose(fid) ;

