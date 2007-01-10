function ok = fmri_svpar(par,parfile,label,trun)
%
% ok = fmri_svpar(par,parfile,<label>,<trun>)
%
% Saves time/stimid to  pardigm file.  If label is included,
% it must be a string matrix with a condition name stored
% on each row.  This name is looked up and appended to each
% row of the paradigm file. trun is the duration of the run.
%
%


%
% fmri_svpar.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

ok = 0;

if(nargin < 2 | nargin > 4)
  msg = 'USAGE: ok = fmri_svpar(par,parfile,<label>,<trun>)';
  qoe(msg);error(msg);
end

nTP = size(par,1);

fid=fopen(deblank(parfile),'w');
if( fid == -1 )
  msg = sprintf('Could not open dof file %s\n',parfile);
  qoe(msg);  error(msg);
end


for n = 1:nTP
  t = par(n,1);
  c = par(n,2);
  if(n < nTP) 
    dt = par(n+1,1) - par(n,1);
  else 
    if(nargin == 4) dt = trun - par(n,1);
    else            dt = -11111;
    end
  end
  if(nargin == 2)
    fprintf(fid,'%6.2f   %2d    %7.4f\n',t,c,dt);
  else
    if(c+1 <= size(label,1))
      fprintf(fid,'%7.3f   %2d     %7.4f   %s\n',t,c,dt,label(c+1,:));
    else
      fprintf(fid,'%7.3f   %2d     %7.4f   %s\n',t,c,dt,'unknown');
    end
  end
end

fclose(fid);
ok = 1;
return;
