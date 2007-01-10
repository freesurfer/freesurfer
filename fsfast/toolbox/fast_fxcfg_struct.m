function fxcfg = fast_fxcfg_struct
% Create an FS-FAST Effects Configuration structure with all 
% the fields. 


%
% fast_fxcfg_struct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.4 $
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

fxcfg.fxtype = '';  % fixed or random
fxcfg.label = '';   % effect name
fxcfg.model = '';   % model string
fxcfg.params = [];  % non-string parameters
fxcfg.sparams = []; % string parameters
fxcfg.npmlist = []; % List of non-parametric matrices

% Cell array of regressor indices. For fixedfx, there is only
% one cell, so regind = fxcfg.regind{1}. For randomfx, the
% cell index indicates the run to use.
fxcfg.regind = []; 

% fxcfg.npmlist(nthrun).M 












