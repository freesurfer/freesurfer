function fxcfg = fast_fxcfg_struct
% Create an FS-FAST Effects Configuration structure with all 
% the fields. 

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












