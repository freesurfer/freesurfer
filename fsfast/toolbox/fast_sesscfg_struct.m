function sesscfg = fast_sesscfg_struct
% Create an FS-FAST Session Configuration structure with all 
% the fields. 

sesscfg.path = '';          % Absolute session path
sesscfg.ntp = [];           % Ntp for each run
sesscfg.runweight = [];     % List of weights for each run
sesscfg.fstemlist = '';     % List of abs path to func stems for each run
sesscfg.evschfilelist = ''; % List of abs path to EvSch Files
sesscfg.evschlist = [];     % Array of event schedules, one for each run

% evschlist(nthrun).evsch
% evsch must have at least 2 columns:
% (1) presentation time
% (2) event id (-1 for tpx)
% (3) presentation weight (optional)

return;







