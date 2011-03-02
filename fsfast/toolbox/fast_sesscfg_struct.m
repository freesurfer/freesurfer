function sesscfg = fast_sesscfg_struct
% Create an FS-FAST Session Configuration structure with all 
% the fields. 


%
% fast_sesscfg_struct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.6 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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







