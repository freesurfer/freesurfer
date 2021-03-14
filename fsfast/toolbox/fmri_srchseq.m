% fmri_srchseq.m
%
% Search for an optimum stimulus sequence.
%
%


%
% fmri_srchseq.m
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

if(0)
  clear all;
  %%%% Optimization Criterion
  % cbe - optimize counter-balancing
  % tr  - minimize expected error in hemodyn response
  Metric = 'cbe';  

  TR = 2.0;

  %% Time for hemodynaic to decay to zero
  TimeWindow = 20.0; 

  % Number of presentations for each condition per run
  nPerCond = [50 75]; % 50 for condition 0, 75 for cond 1

  % Number of paradigm files to create (ie, number of runs)
  % Create many more than you need as this does not slow
  % the search down by much.
  nRuns = 10; 

  %% Stem of the name of the output paradigm files
  %% Full name is parstemXXX.dat where XXX goes from
  %% 000 to nRuns-1
  parstem = 'par';

  %% Number of sequences overwhich to search for the
  %% optimum.  This should be 10^5 to 10^7.  This may
  %% need to be run over night.
  nSearch = 100;

  %% Seed for random number generator %%%
  RSeed = 53;

  nSwap = 30;
end

%%%%%%%%%-------------------------------%%%%%%%%%%%%%%%
% No changes need to be made below this line
%%%%%%%%%-------------------------------%%%%%%%%%%%%%%%
fprintf('\n$Id: fmri_srchseq.m,v 1.3 2011/03/02 00:04:06 nicks Exp $\n');

RunDuration = TR*Ntp;
nStim = floor(RunDuration/TER);
n0 = nStim - sum(nPerCond);
if(n0 < 0)
  msg = 'ERROR: too many stimuli for given run duration';
  qoe(msg);error(msg);
end
nPerCond = [n0 nPerCond];

if(nSwap < 1) nSwap = sum(nPerCond); end

if(nSwap > sum(nPerCond))
  msg = sprintf('\nERROR: nSwap (%d) > sum(nPerCond) (%d)',nSwap,sum(nPerCond));
  qoe(msg);error(msg);
end

Nh = round(TimeWindow/TER);
if(Nh < 1)
  msg = sprintf('TimeWindow (%g) must be >= TER (%g)',TimeWindow,TER);
  qoe(msg);error(msg);
end

rand('state',RSeed);
nCond = length(nPerCond);  % Number of Conditions
nNNCond = nCond - 1;       % Number of Non-Null Conditions

nKeep = nRuns;
timewindow = TimeWindow;
Nc = nCond;
Nnnc = nNNCond;
Ntot = Nh*Nnnc;

nSeq = sum(nPerCond);

if( strcmp(Metric,'tr') | strcmp(Metric,'minerr') )
  if(Ntp <= Ntot)
    fprintf('ERROR: Too many estimates.');
    fprintf('With %d conditions and %d time points, Nh = %d max\n',...
        Nc,Ntp,floor(Ntp/Nc));
    fprintf('Reduce timewindow from %g to %g\n',TimeWindow,TER*floor(Ntp/Nc));
    msg = '';
    qoe(msg);error(msg);
  end
end

[idealcbpm idealxtx] = fmri_idealcbpm([nPerCond],Nh);
tr_ideal = trace(inv(idealxtx));

seqbase = fmri_seqbase(nPerCond);
[Mss Rss]   = fmri_subsampmat(TR,Ntp,TER);

fprintf(1,'Nc   = %d\n',Nc);
fprintf(1,'nSeq = %d\n',nSeq);
fprintf(1,'Nh   = %d\n',Nh);
fprintf(1,'Rss  = %d\n',Rss);
fprintf(1,'Ntot = %d\n',Ntot);

fprintf(1,'size xtx = %d %d\n',size(idealxtx));

fprintf(1,'Starting Search over %d items\n',nSearch);

if( strcmp(Metric,'tr') | strcmp(Metric,'minerr') )
  tic;
  if(Rss == 1) 
    [seq tr travg trstd] = fmri_minseqtr(seqbase,Nh,nSearch,nKeep,nSwap);
  else
    [seq tr travg trstd] = fmri_minseqtr(seqbase,Nh,nSearch,nKeep,nSwap,Mss);
  end
  tSrch = toc;
  fprintf(1,'Search Time = %g\n',toc);
elseif (strcmp(Metric,'cbe') )
  tic;
  [seq cbe cbeavg cbestd] = ...
       fmri_minseqcbe(seqbase,Nh,nSearch,nKeep,idealxtx,nSwap);
  tSrch = toc;
  fprintf(1,'Search Time = %g\n',toc);
else
  msg = sprintf('Unrecognized metric %s, use cbe or tr');
  qoe(msg);error(msg);
end

scm  = fmri_seq2scm(seq,Nh);
for k = 1:nKeep,
  x = scm(:,:,k);
  xtx = x'*Mss'*Mss*x; %'
  cbe(k) = mean(abs(reshape1d(idealxtx-xtx)));
  tr(k)  = trace(inv(xtx));
end

t = TER * [0:nSeq-1]';
l = reshape1d(([1 1 1]' * [1:8]));

for k = 1:nKeep,
  par = [t seq(:,k)];

  parfile = sprintf('%s-%03d.dat',parstem,k);
  fmri_svpar(par,parfile);

end

%%%%%%%%%%% --------- save the info file ------- %%%%%
infofile = sprintf('%s-srch.info',parstem);
fid = fopen(infofile,'w');
if(fid == -1)
  msg = sprintf('Could not open %s for writing',infofile);
  qoe(msg);error(msg);
end

fprintf(fid,'Metric %s\n',Metric);
fprintf(fid,'nSearch %d\n',nSearch);
fprintf(fid,'RSeed   %d\n',RSeed);
fprintf(fid,'tSearch %g sec = %g hrs\n',tSrch,tSrch/3600);
fprintf(fid,'nSwap %d\n',nSwap);
fprintf(fid,'Nh %d\n',Nh);
fprintf(fid,'Ntp %d\n',Ntp);
fprintf(fid,'TR %g\n',TR);
fprintf(fid,'TER %g\n',TER);
fprintf(fid,'nPerCond ');
fprintf(fid,'%3d ',nPerCond);
fprintf(fid,'\n');
fprintf(fid,'Trace of Min CBError %g\n',tr_ideal);
fprintf(fid,'Run   Trace    CBError\n');
for k = 1:nKeep,
  fprintf(fid,'%3d %8.4f %8.4f\n',k,tr(k),cbe(k));
end
fclose(fid);

return;
