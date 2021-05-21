function qoe(msg)
%
% qoe(<msg>)
%
% When there is a global variable defined called QuitOnError
% and it is 1, then matlab is quit after printing out the
% backtrace.  To use, declare "global QuitOnError" and set
% to either 0 or 1.
%
% The main use of this function is when a matlab script is being 
% executed non-interactively.  If the error() function were called
% instead, then matlab would hand control back to a non-existent
% terminal causing the script to crash very ungracefully.  Remember
% declare global QuitOnError when starting up.
%
%


%
% qoe.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

global QuitOnError

%fprintf('qoe: QuitOnError = %d\n',QuitOnError);

if(exist('QuitOnError'))
  if(~isempty(QuitOnError))
    if(QuitOnError==1)
      %st = dbstack; % get the call stack 
      st = []; % get the call stack 
      nst = size(st,1);
      for n = 2:nst,
        fprintf(2,'%s Line %d\n',...
          getfield(st(n,:),'name'),getfield(st(n,:),'line'));
      end
      if(nargin==1) fprintf(2,'%s\n',msg); end
      fprintf('quiting matlab\n');
      quit force;
    end
  end
end
return;

