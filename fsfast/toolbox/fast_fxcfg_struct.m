function fxcfg = fast_fxcfg_struct
% Create an FS-FAST Effects Configuration structure with all 
% the fields. 

fxcfg.fxtype = '';  % fixed or random
fxcfg.label = '';   % effect name
fxcfg.model = '';   % model string
fxcfg.params = [];  % non-string parameters
fxcfg.sparams = []; % string parameters
fxcfg.npmlist = []; % List of non-parametric matrices

% fxcfg.npmlist(nthrun).M 












