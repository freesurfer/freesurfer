function sesscfg = fast_sesscfg_struct
% Create an FS-FAST Session Configuration structure with all 
% the fields. 


%
% fast_sesscfg_struct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

sesscfg.ntp = [];           % Ntp for each run
sesscfg.runlist = [];
sesscfg.runweight = [];     % List of weights for each run
sesscfg.fstemlist = '';     % List of abs path to func stems for each run
sesscfg.evschfilelist = ''; % List of abs path to EvSch Files (needed?)
sesscfg.evschlist = [];     % Array of event schedules, one for each run
sesscfg.volsize = [];       % rows, cols, slices

% evschlist(nthrun).evsch
% evsch must have at least 2 columns:
% (1) presentation time
% (2) event id (-1 for tpx)
% (3) presentation weight (optional)

return;







