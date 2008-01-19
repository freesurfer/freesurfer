function ev = flac_ev_parse(tline)
% ev = fast_ev_parse(<tline>)
%
% Parses an EV line from a FLAC file. If a tline is
% not provided, then returns an empty ev struct.
%
% EV EVName ModelName Type <parameters>
%
%


%
% flac_ev_parse.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2008/01/19 22:32:05 $
%    $Revision: 1.15 $
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
ev.segstem    = '';  % used for selfregseg
ev.varreg     = 0;   % EV has a variable number of regressors
% VarRegEVs must be nuissance!
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
 case {'gamma'} % Gamma HRF regressor
  % 5 parameters: delay dispersion alpha nderiv dpsd
  % EV Probe gamma task 2.25 1.25 2 0 .1
  [ev.stf c] = sscanfitem(tline,5);
  if(c ~= 1) 
    fprintf('Format error: %s: stf\n',ev.model); 
    fprintf('Expected something like: EV Probe gamma task  2.25 1.25 2  0 .1\n');
    ev=[]; 
    return; 
  end

  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error: %s: delay\n',ev.model); ev=[]; return; end
  ev.params(1) = sscanf(item,'%f',1); % delay (sec)

  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error: %s: dispersion\n',ev.model); ev=[]; return; end
  ev.params(2) = sscanf(item,'%f',1); % dispersion (sec)

  [item c] = sscanfitem(tline,8);
  if(c ~= 1) fprintf('Format error: %s: alpha\n',ev.model); ev=[]; return; end
  ev.params(3) = sscanf(item,'%f',1); % alpha

  [item c] = sscanfitem(tline,9);
  if(c ~= 1) fprintf('Format error: %s: nderiv\n',ev.model); ev=[]; return; end
  ev.params(4) = sscanf(item,'%f',1); % nderiv

  [item c] = sscanfitem(tline,10);
  if(c ~= 1) fprintf('Format error: %s: dpsd\n',ev.model); ev=[]; return; end
  ev.params(5) = sscanf(item,'%f',1); % dpsd (sec)

  ev.psdwin = [0 40 ev.params(5)];
  ev.npsd = round((ev.psdwin(2)-ev.psdwin(1))/ev.psdwin(3));
  ev.ishrf = 1;  
  ev.nreg = 1+ev.params(4); % 1+nderiv
  
 %--------------------------------------------
 case {'fslgamma'} % FSL Gamma HRF regressor
  % 5 parameters: phase std meanlag nderiv dpsd
  % EV Probe fslgamma task stfile 0 3 6 0 .1
  [ev.stf c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error: %s: stf\n',ev.model); ev=[]; return; end

  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error: %s: phase\n',ev.model); ev=[]; return; end
  ev.params(1) = sscanf(item,'%f',1); % phase (sec)

  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error: %s: stddev\n',ev.model); ev=[]; return; end
  ev.params(2) = sscanf(item,'%f',1); % stddev (sec)

  [item c] = sscanfitem(tline,8);
  if(c ~= 1) fprintf('Format error: %s: meanlag\n',ev.model); ev=[]; return; end
  ev.params(3) = sscanf(item,'%f',1); % meanlag

  [item c] = sscanfitem(tline,9);
  if(c ~= 1) fprintf('Format error: %s: nderiv\n',ev.model); ev=[]; return; end
  ev.params(4) = sscanf(item,'%f',1); % nderiv

  [item c] = sscanfitem(tline,10);
  if(c ~= 1) fprintf('Format error: %s: dpsd\n',ev.model); ev=[]; return; end
  ev.params(5) = sscanf(item,'%f',1); % dpsd (sec)

  ev.psdwin = [0 40 ev.params(5)];
  ev.npsd = round((ev.psdwin(2)-ev.psdwin(1))/ev.psdwin(3));
  ev.ishrf = 1;  
  ev.nreg = 1+ev.params(4); % 1+nderiv
  
 %--------------------------------------------
 case {'selfregseg'} % Segmentation-based Self-Regressor
   % Variable number of params
   % EV CSF selfreg nuis segstem ResidFlag nPCA Id1 ... IdN
   % EV CSF selfreg nuis aparc+aseg 0 10   4 5 43 44 14 15 72
   % nregressors = nPCA
   [item c] = sscanfitem(tline,5);
   if(c ~= 1) fprintf('Format error: %s: SegStem\n',ev.model); ev=[]; return; end
   ev.segstem = item; % segmentation volume
   
   [item c] = sscanfitem(tline,6);
   if(c ~= 1) fprintf('Format error: %s: ResidFlag\n',ev.model); ev=[]; return; end
   ev.params(1) = sscanf(item,'%d',1); % resid flag
   
   [item c] = sscanfitem(tline,7);
   if(c ~= 1) fprintf('Format error: %s: nPCA\n',ev.model); ev=[]; return; end
   ev.params(2) = sscanf(item,'%d',1); % nPCA

   nthid = 0;
   while(1)
     [item c] = sscanfitem(tline,8+nthid);
     if(c ~= 1) break; end
     ev.params(3+nthid) = sscanf(item,'%d',1); % SegId
     nthid = nthid + 1;
   end
   if(nthid == 0) fprintf('Format error: %s: SegId\n',ev.model); ev=[]; return; end
   ev.nreg = ev.params(2);
   ev.ishrf = 0;
   
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
  % EV MCReg nonpar nuis mcextreg 3  
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
 case {'asl'} % Just ones on the odds or the evens
  % 1 parameters
  % EV ASL asl task odd-or-even
  [item c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error: %s: odd-or-even\n',ev.model); ev=[]; return; end
  if(strcmp(item,'tag-first'))   ev.params(1) = 0; end
  if(strcmp(item,'tag-second'))  ev.params(1) = 1; end
  ev.nreg = 1;
  ev.ishrf = 0;  

 %--------------------------------------------
 case {'aareg'} % 
  % 3 parameters
  % EV AAREG aareg nuis nkeep fmax fdelta
  [item c] = sscanfitem(tline,5);
  if(c ~= 1) fprintf('Format error: %s\n',tline); ev=[]; return; end
  ev.params(1) = sscanf(item,'%d',1); % nkeep
  ev.nreg = ev.params(1);

  [item c] = sscanfitem(tline,6);
  if(c ~= 1) fprintf('Format error: %s\n',tline); ev=[]; return; end
  ev.params(2) = sscanf(item,'%f',1); % fmax

  [item c] = sscanfitem(tline,7);
  if(c ~= 1) fprintf('Format error: %s\n',tline); ev=[]; return; end
  ev.params(3) = sscanf(item,'%f',1); % fdelta
  
  ev.ishrf = 0;  

 %--------------------------------------------  
 otherwise
  fprintf('ERROR: flac_ev_parse: model %s unrecoginized\n',ev.model);
  ev = [];
  return;
  
end % switch


return;
