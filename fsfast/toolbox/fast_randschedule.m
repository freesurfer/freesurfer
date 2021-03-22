function [par, MeanISI] = fast_randschedule(Npc,Tpc,RunDuration,TPreScan,Nruns,Tres)
% [par MeanISI] = fast_randschedule(Npc,Tpc,RunDuration,TPreScan,Nruns,Tres)
%
% Npc = vector of number of presentations of each stimulus type
% Tpc = duration of the presentation of each stimulus type
% RunDuration = time of the run (Ntp*TR+PreScan)
% TPreScan = time during which stimuli should be presented before
%   any images are stored
% Nruns = number of runs to create
% Tres = temporal resolution (-1 for infinite)
%
% par = Ns X 2 X Nruns, where Ns is the number of stimulus presentations
%   (ie, Ns = sum(Npc)), 2 = [onset-time id]


%
% fast_randschedule.m
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
MeanISI = [];

if(nargin ~= 6)
  fprintf('USAGE: par = fast_randschedule(Npc,Tpc,RunDuration,TPreScan,Nruns,Tres)\n');
  return;
end

Nc = length(Npc);
if(Nc ~= length(Tpc))
  fprintf('ERROR: length of Npc does not equal length of Tpc\n');
  return;
end

Nstimtot = sum(Npc);
Ttot = RunDuration + TPreScan;
Tstimtot = sum(Npc.*Tpc);
if(Ttot <= Tstimtot)
  fprintf('ERROR: Total scan time = %g < Total stimulation time %g\n',...
          Ttot,Tstimtot);
  return;
end

MeanISI = Ttot/(Nstimtot-1);

Tnulltot = Ttot - Tstimtot;
Tnullavg = Tnulltot/Nstimtot;

StimIdBase = [];
for n = 1:Nc
  StimIdBase = [StimIdBase; n*ones(Npc(n),1)];
end

par = zeros(Nstimtot,2,Nruns);

for run = 1:Nruns

  Tnull = rande([Nstimtot 1]);
  Tnull = Tnulltot*Tnull/sum(Tnull);
  StimId = StimIdBase(randperm(Nstimtot));
  StimDur = reshape1d(Tpc(StimId));
  StimOnset = cumsum(StimDur + Tnull);
  StimOnset = [0; StimOnset(1:Nstimtot-1)];
  StimOnset = StimOnset - TPreScan;

  if(Tres > 0)
    StimOnset = Tres * round(StimOnset/Tres); 
  end

  par(:,1,run) = StimOnset;
  par(:,2,run) = StimId;

end


return;
