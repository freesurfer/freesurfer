function sesscfg = fast_sesscfg_struct
% Create an FS-FAST Session Configuration structure with all 
% the fields. 

sesscfg.path = '';          % Absolute session path
sesscfg.ntp = [];           % Ntp for each run
sesscfg.runweight = [];     % List of weights for each run
sesscfg.fstemlist = '';     % List of abs path to func stems for each run
sesscfg.evschfilelist = ''; % List of abs path to EvSch Files
sesscfg.evsch = [];         % Array of event schedules, one for each run

% evsch.tpres(nthrun)
% evsch.evidpres(nthrun)
% evsch.wpres(nthrun)

return;







