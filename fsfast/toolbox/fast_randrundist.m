function [Npcrun, TimeOut, NTries] = ...
   fast_randrundist(Npc,Tpc,RunDuration,Nruns,NTriesMax)
% [Npcrun, TimeOut, NTries] = 
%    fast_randrundist(Npc,Tpc,RunDuration,Nruns,NTriesMax)
%
% Distributes number of stimulus presentations across runs
% according to a semi-random, iterative formula. 
%
% Npc = Number of presentations per condition across all runs, 
%       excluding null stimulus. Size must match Tpc.
% Tpc = Duration of each presentation of each condition.
% NTriesMax = maximum number of iterations
%
% Npcrun = (Nc-by-Nruns) matrix of number of stimuli of each 
%          type to present on each run.
% TimeOut = 1 for a timeout error (Npcrun will also be []).
% NTries  = number of iterations.


%
% fast_randrundist.m
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

Npcrun = [];
TimeOut = 0;
NTries = -1;

if(nargin ~= 5)
  fprintf('[Npcrun TimeOut NTries] = fast_randrundist(Npc,Tpc,RunDuration,Nruns,NTriesMax)\n');
  return;
end

Npc = reshape1d(Npc);
Tpc = reshape1d(Tpc);
Nc = length(Npc);

if(length(Tpc) ~= Nc)
  fprintf('ERROR: size of Npc does not match size of Tpc\n');
 return;
end


TTot = Nruns*RunDuration;
TStimTot = sum(Npc.*Tpc);
if(TStimTot >= TTot)
 fprintf('ERROR: total stimulation time %g is >= than total time %g\n',...
         TStimTot, TTot);
 return;
end

NpcRepMat = repmat(Npc,[1 Nruns]);
TpcRepMat = repmat(Tpc,[1 Nruns]);

stop = 0;
NTries = 0;
while(~stop)

  r = rand(Nc,Nruns);
  rsum = sum(r,2);

  r = r./repmat(rsum,[1 Nruns]);

  Npcrun = round( r .* NpcRepMat);

  TStimRun = sum(Npcrun .* TpcRepMat);

  nover = length(find(TStimRun >= RunDuration));

  NTries = NTries + 1;

  if(nover == 0)
    stop = 1;
  elseif(NTries > NTriesMax)
    stop = 1;
    TimeOut = 1;
    Npcrun = [];
  end

end


return;
