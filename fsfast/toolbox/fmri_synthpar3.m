function [par,t0avg] = fmri_synthpar3(npercond,tpercond,nruns,trun,tres,tprescan)
% [par,t0avg] = fmri_synthpar3(npercond,tpercond,nruns,trun,tres,tprescan)
%
% npercond - number of presentations per NON-NULL condition, all runs
% tpercond - stimulus duration of each NON-NULL condition
% nruns    - number of runs 
% trun     - duration of each run
% tres     - temporal resolution (0 for infinite)
% tprescan - time of first stimulus prior to scanning (not included
%            in trun).
%


%
% fmri_synthpar3.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

par = [];
t0avg = 0;

if(nargin ~= 6)
  msg = 'USAGE: par = fmri_synthpar3(npercond,tpercond,nruns,trun,tres,tprescan)';
  qoe(msg);error(msg);
end

Nconds = length(npercond);
if(length(tpercond) ~= Nconds)
  msg = sprintf('Dimension of npercond (%d) does not equal tpercond (%d)',...
                Nconds, length(tpercond));
  qoe(msg);error(msg);
end

if(tres > 0)
  tpercond = tres*ceil(tpercond/tres);
end

nstim    = sum(npercond);              % total number of stimuli in session
truntot  = trun + tprescan;            % duration of each run
tsession = nruns*truntot;              % total session time
tstimtot = sum(tpercond .* npercond);  % total stimulus time in session
if(tstimtot > tsession)
  msg = sprintf('tstimtot (%g) > tession (%g)',tstimtot,tsession);
  qoe(msg);error(msg);
end

t0tot = tsession - tstimtot; % total non-stimulus time (ie, null)
t0avg = t0tot/(nstim-1);     % average null between each event

%% Create null events (variable duration) %%
if(tres <= 0)
  p0 = rand(nstim,1);    % create random vector, one elem for each pres
  p0 = p0/sum(p0);       % normalize 
  t0 = t0tot*p0;         % vector of null times between each event
else
  n0 = nulldist(t0tot,nstim,tres);
  t0 = tres*n0;
end

%%% Create a list of non-null events and their durations %%%%
eventlist = [];
eventtimelist = [];
for cond = 1:Nconds,
  eventlist = [eventlist; cond*ones(npercond(cond),1)];
  eventtimelist = [eventtimelist; tpercond(cond)*ones(npercond(cond),1)];
end

%%% Randomize list of non-null events (same rand for their durations) %%%%
rp = randperm(nstim);
eventlist = eventlist(rp);
eventtimelist = eventtimelist(rp);

%tt = eventtimelist + t0;
%ctt = [cumsum([-tprescan tt])];

%%% Go through each run %%%
nevent = 1;
for run = 1:nruns

  %%% Go through time on each run %%%
  t = -tprescan;
  runevent = 1;
  sumt0 = 0;
  nt00 = 0;
  while(t < trun)
    event = eventlist(nevent);   % Get the event code, 
                                 % nevent --> number of non-null session events

    par(runevent,1,run) = t;     % Set the event onset time in the paradigm
    par(runevent,2,run) = event; % Set the event code in the paradigm

    t = t + tpercond(event);     % This is the end of the event
    runevent = runevent + 1;     % Event number inside the run (incl null)

    sumt0 = sumt0 + t0(nevent);

    %% Insert the null event %%
    if(t0(nevent) > 0)
      par(runevent,1,run) = t; % Set the null-event onset time in the paradigm
      par(runevent,2,run) = 0; % Code is always 0 
      t = t + t0(nevent);      % This is the end of the null-event
      runevent = runevent + 1; % Event number inside the run (incl null)
    else
      nt00 = nt00 + 1;
    end

    nevent = nevent + 1; % update the number of sesssion non-null events so far


    %% Check to see if all the session events have been used %%
    if(nevent > nstim) 
      %% Add a null event at the end of the run %%
      %% This is a hack for when paradigms are a different length %%
      par(runevent,1,run) = t; 
      par(runevent,2,run) = 0; 
      tresidual = (t-trun);
      %fprintf('run = %2d, tresidual = %g\n',run,tresidual);
      return; 
    end

  end % while t < trun loop

  % Check whether the stimulus extended beyond the end of the run. %
  % If so, account for the time by adding to future null events
  tresidual = t-trun;
  %fprintf('run = %2d, tresidual = %g\n',run,tresidual);
  n0residual = floor(tresidual/tres);
  if(n0residual > 0)
    dn = min(n0residual,nstim-nevent);
    % rp = nevent + randperm(nstim-nevent);
    rp = nevent + randperm(dn);
    t0(rp(1:dn)) = t0(rp(1:dn)) + tres;
  end


  %% Add a null event at the end of the run %%
  %% This is a hack for when paradigms are a different length %%
  par(runevent,1,run) = t; 
  par(runevent,2,run) = 0; 

end % loop over runs

if(nevent < nstim) 
  fprintf('WARNING: fmri_synthpar3: nused = %d, nstim = %d\n',...
          nevent,nstim) ;
end


return;

%-------------------------------------------------------%
function n0 = nulldist(t0tot,nstim,tres)
%
% n0 = nulldist(t0tot,nstim,tres)
%

n0tot = floor(t0tot/tres); % number of minimum-duration null events

lmax = 20;
lm = 0;
dn = 1;
while(dn ~= 0)

  %randn('state')
  p0 = rande(nstim);     % create random (exp dist) vector, one elem for each pres
  p0 = p0/sum(p0);       % normalize 
  t0 = t0tot*p0;         % vector of null times between each event

  n0 = round(n0tot*p0);  % round to number of minimum-duration events (MDEs)
  dn = n0tot-sum(n0);    % sum(n0) = actual number of MDEs

  if(dn > 0)
    % there are not enough actual MDEs, add the difference to random
    % slots in n0.
    rp = randperm(nstim);
    n0(rp(1:dn)) = n0(rp(1:dn)) + 1;
  end

  dn = n0tot-sum(n0);

  lm = lm + 1;
  if(lm > lmax)
    msg = sprintf('fmri_synthpar3: timed out');
    qoe(msg);error(msg);
  end

end

return;
%----------------------------------------------------------%
