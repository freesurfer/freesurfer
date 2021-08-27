function RM = fmri_mrestriction(TestType,nH, nC, AC, CC, nHTest)
%
% RM = fmri_mrestriction(TestType, nH, nC, AC, CC, <nHTest>)
%
% Creates a restriction matrix 
%
% TestType - String (case insensitive) indicating the test to perform:
%   t  - t-Test, RM will be a row vector
%   tm - t-Test, different p for each delay point (RM same as Fc)
%   Fm - F-Test, different p for each delay point
%   F0 - RM same as t-Test
%   Fd - F-test, different row for each delay, all conditions
%        for each delay on the same row
%   Fc - F-test, different row for each condition, all delays
%        for a condition on the same row
%   Fcd - F-test, different row for each condition for each delay
% nH - total number of elements in HDIR
% nC - total number of stimulus conditions (including fixation)
% AC - list of active conditions
% CC - list of control conditions
% nHTest - list of HDIR components to test (default: all).
%
% Eg: fmri_mrestriction('t',10, 6, [1 3], [2 5], [3 6:9])
% Generates an RM for testing conditions (1+3)-(2+5) using
% components 3,6,7,8, and 9 in the HDIR.  Conditions 0 and 4 
% are not tested. There are 10 components in the HDIR, and 5 
% stimulus conditions.
%
%
%


%
% fmri_mrestriction.m
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
if(nargin ~= 5 & nargin ~= 6)
  msg = 'USAGE: RM = fmri_mrestriction(TestType, nH, nC, AC, CC, <nHTest>';
  qoe(msg);error(msg);
end

% Check the Test Type %
if( isempty( strmatch(upper(TestType),{'T ','TM','FM','F0','FD','FC','FCD','FDC'},'exact')))
  msg = sprintf('Unkown TestType %s',TestType);
  qoe(msg);error(msg);
end  

%% check that the input parameters are correct %%
if(nH < 1)
  msg = 'nH must be greater than 0';
  qoe(msg);error(msg);
end
if(nC < 1)
  msg = 'Number of conditions must be greater than 0';
  qoe(msg);error(msg);
end
if(length(find(AC<0 | AC > nC)))
  msg = 'Invalid condition number in AC';
  qoe(msg);error(msg);
end
if(length(find(CC<0 | CC > nC)))
  msg = 'Invalid condition number in CC';
  qoe(msg);error(msg);
end

% Test that the same condition is not in both AC and CC %
C = [AC CC];
for n = 1:nC,
  if(length(find(C==n)) > 1)
    msg = 'Same condition is in both AC and CC';
    qoe(msg);error(msg);
  end
end

% Set defaults for nHTest, if needed, else test range %
if(nargin == 5) nHTest = [1:nH];
else
  if(length(find(nHTest<1 | nHTest > nH)))
    msg = 'Invalid nHTest Component';
    qoe(msg);error(msg);
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
switch(upper(TestType))
  case {'T','F0'}, 
    RM = fmri_vrestriction(nH,nC,AC,CC,nHTest);
  case {'FD','TM','FM'},     
    for n = 1:length(nHTest),
      RV = fmri_vrestriction(nH,nC,AC,CC,nHTest(n));
      RM(n,:) = reshape(RV', [1 prod(size(RV))]);
    end
  case {'FC'},
    for n = 1:length(iACnz),
      RV = fmri_vrestriction(nH,nC,AC(n),0,nHTest);
      RM(n,:) = reshape(RV', [1 prod(size(RV))]);
    end
    for n = 1:length(iCCnz),
      RV = fmri_vrestriction(nH,nC,0,CC(n),nHTest);
      RM(n+length(AC),:) = reshape(RV', [1 prod(size(RV))]);
    end
  case {'FCD','FDC'},
    RM = [];
    for n = 1:length(iACnz),
      RV = fmri_mrestriction('Fd',nH,nC,AC(n),0,nHTest);
      RM = [RM; RV;];
    end
    for n = 1:length(iCCnz),
      RV = fmri_mrestriction('Fd',nH,nC,0,CC(n),nHTest);
      RM = [RM; RV;];
    end
end

return;

