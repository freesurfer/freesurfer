function ev = flac_ev_parse(tline)
% ev = fast_ev_parse(<tline>)
%
% Parses an EV line from a FLAC file. If a tline is
% not provided, then returns an empty ev struct.
%
% EV EVName ModelName Type <parameters>
%
% $Id: flac_ev_parse.m,v 1.6 2004/11/01 04:38:08 greve Exp $

ev = [];
if(nargin > 1)
  fprintf('ev = fast_ev_parse(<tline>)\n');
  return;
end

ev.name       = '';  % Name of the EV
ev.model      = [];  % Name of the model
ev.type       = [];  % task or nuis
ev.params     = [];  % params
ev.ishrf      = 0;   % Flag for HRFs
ev.psdwin     = [];  % HRF psd window [min d max]
ev.npsd       = [];  % Number of items in psd win
ev.stf        = '';  % Stimulus timing file name
ev.st         = [];  % Stimulus timing [tonset duration weight]
ev.Xfir       = [];  % FIR matrix for HRFs
ev.Xirf       = [];  % IRF matrix for HRFs
ev.nonparname = '';  % Name of non-parametric regressor
ev.X          = [];  % Design matrix for this EV
ev.nreg       = [];  % number of regressors
if(nargin == 0)  return; end

% Read in the name
[ev.name c] = sscanf(tline,'%*s %s',1);
if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  
% Read in the model type
[ev.model c] = sscanf(tline,'%*s %*s %s',1);
if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end

% Read in the EV type (should be task or nuis)
[ev.type c] = sscanf(tline,'%*s %*s %*s %s',1);
if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
if(~strcmp(ev.type,'task') & ~strcmp(ev.type,'nuis'))
  fprintf('ERROR: %s,%s: EVType is %s, must be nuis or task\n',...
	  ev.name, ev.model, ev.type);
  ev = [];
  return;
end

switch (ev.model)
 
 %--------------------------------------------
 case {'baseline'} % Baseline/Mean offset
  % 0 parameters:
  % EV  Baseline baseline nuiss
  ev.nreg = 1;
  ev.ishrf = 0;  
  
 %--------------------------------------------
 case {'polynomial'} % Polynomial
  % 1 parameter: order
  % EV   Poly   polynomial nuiss 2
  % Note: does not include mean offset
  [item c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.params(1) = sscanf(item,'%d',1); % order
  ev.nreg = ev.params(1);
  ev.ishrf = 0;  
  
 %--------------------------------------------
 case {'fir'} % FIR
  % 4 parameters: stf, psdmin, dpsd, psdmax 
  % EV SMFIR  fir  task  sm.stf -4 20 .5
  [ev.stf c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  
  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.psdwin(1) = sscanf(item,'%f',1); % psdmin
  
  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.psdwin(2) = sscanf(item,'%f',1); % psdmax
  
  [item c] = sscanfitem(tline,8);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.psdwin(3) = sscanf(item,'%f',1); % dpsd

  ev.npsd = round((ev.psdwin(2)-ev.psdwin(1))/ev.psdwin(3));
  ev.nreg = ev.npsd;
  ev.ishrf = 1;  

 %--------------------------------------------
 case {'spmhrf'} % SPM HRF regressor
   % 3 parameters: stf, nderiv, dpsd
  [ev.stf c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  
  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.params(1) = sscanf(item,'%d',1); % nderiv
  
  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.params(2) = sscanf(item,'%f',1); % dpsd
  
  ev.psdwin = [0 40 ev.params(2)];
  ev.npsd = round((ev.psdwin(2)-ev.psdwin(1))/ev.psdwin(3));
  ev.ishrf = 1;    
  ev.nreg = ev.params(1) + 1;
  

 %--------------------------------------------
 case {'fourier'} % Fourier regressor
  % 3 parameters: period nharmonics tdelay
  % EV SMPer fourier task 30 2 3
  [item c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error: %s: period\n',ev.model); ev=[]; return; end
  ev.params(1) = sscanf(item,'%f',1); % period (sec)

  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error: %s: nharm\n',ev.model); ev=[]; return; end
  ev.params(2) = sscanf(item,'%d',1); % nharmonics

  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error: %s: tdelay\n',ev.model); ev=[]; return; end
  ev.params(3) = sscanf(item,'%f',1); % tdelay (sec)

  ev.nreg = 2*(ev.params(2)+1); % 2*(nharm+1)
  ev.ishrf = 0;  
  
 %--------------------------------------------
 case {'nonpar'} % Non-parametric regressor
  % 2 parameters: name, ncols
  % EV MCReg nonpar nuiss mcextreg 3  
  [ev.nonparname c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end

  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error\n'); ev=[]; return; end
  ev.params(1) = sscanf(item,'%d',1); % ncols
  ev.nreg = ev.params(1);
  ev.ishrf = 0;  

 %--------------------------------------------
 case {'nyquist'} % Temporal nyquist
  % 0 parameters
  % EV Nyq nyquist nuiss 
  ev.nreg = 1;
  ev.ishrf = 0;  

 %--------------------------------------------  
 otherwise
  fprintf('ERROR: model %s unrecoginized\n');
  ev = [];
  return;
  
end % switch


return;
