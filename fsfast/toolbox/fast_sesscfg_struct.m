function sesscfg = fast_sesscfg_struct
% Create an FS-FAST Session Configuration structure with all 
% the fields. 

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







