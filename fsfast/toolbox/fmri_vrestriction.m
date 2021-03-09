function RV = fmri_vrestriction(nH, nC, AC, CC, nHTest)
%
% RV = fmri_vrestriction(nH, nC, AC, CC, nHTest)
%
% Creates a restriction vector of dimension 1 X nC*nH
%
% nH - total number of elements in HDIR
% nC - total number of stimulus conditions (including fixation)
% AC - list of active conditions
% CC - list of control conditions
% nHTest - list of HDIR components to test (default: all).
%
% Eg: fmri_vrestriction(10, 6, [1 3], [2 5], [3 6:9])
% Generates an RV for testing conditions (1+3)-(2+5) using
% components 3,6,7,8, and 9 in the HDIR.  Conditions 0 and 4 
% are not tested. There are 10 components in the HDIR, and 5 
% stimulus conditions.
%
%
%


%
% fmri_vrestriction.m
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

% Check for the correct number of arguments %
if(nargin ~= 4 & nargin ~= 5)
  error('USAGE: fmri_vrestriction(nH, nC, AC, CC, <nHTest>');
  return;
end

%% check that the input parameters are correct %%
if(nH < 1)
  msg = 'nH must be greater than 0';
  qoe(msg); error(msg);
end
if(nC < 1)
  msg = 'Number of conditions must be greater than 0';
  qoe(msg); error(msg);
end
if(length(find(AC<0 | AC > nC-1)))
  msg = 'Invalid condition number in AC';
  qoe(msg); error(msg);
end
if(length(find(CC<0 | CC > nC-1)))
  msg = 'Invalid condition number in CC';
  qoe(msg); error(msg);
end

% Test that the same condition is not in both AC and CC %
C = [AC CC];
for n = 1:nC,
  if(length(find(C==n)) > 1)
    msg = 'Same condition is in both AC and CC';
    qoe(msg); error(msg);
  end
end

% Set defaults for nHTest, if needed, else test range %
if(nargin == 4) nHTest = [1:nH];
else
  if(length(find(nHTest<1 | nHTest > nH)))
    msg = 'Invalid nHTest Component';
    qoe(msg); error(msg);
  end
end

% Strip Condition 0 from AC (if there) %
iACnz = find(AC ~= 0);
if( ~ isempty(iACnz) )
  AC = AC(iACnz);
end
% Strip Condition 0 from CC (if there) %
iCCnz = find(CC ~= 0);
if( ~ isempty(iCCnz) )
  CC = CC(iCCnz);
end

%% -------- Generate the restriction matrix -------- %%
RV = zeros(nC-1,nH); % -1 excludes fixation
if( ~ isempty(iACnz) )
  RV(AC,nHTest) =  ones(length(AC),length(nHTest));
end
if( ~ isempty(iCCnz) )
  RV(CC,nHTest) = -ones(length(CC),length(nHTest));
end
RV = reshape(RV', [1 prod(size(RV))]);

return;

