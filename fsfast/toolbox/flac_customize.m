function flacnew = flac_customize(flac)
% flacnew = flac_customize(flac)
%
% Gets number of time points for the given run, loads all stimulus
% timings, and loads all nonpar matrices. 
%
% Caller must set: flac.sess and flac.nthrun
%
% See flac_desmtx for how the design matrices are built.
%
%


%
% flac_customize.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/06/15 23:02:36 $
%    $Revision: 1.32 $
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

flacnew = [];
if(nargin ~= 1)
 fprintf('flacnew = flac_customize(flac)\n');
 return;
end
flacnew = flac;

% Construct path names
fsdpath = sprintf('%s/%s',flac.sess,flac.fsd);
flacnew.runlist = fast_runlist(fsdpath,flac.runlistfile);
if(isempty(flacnew.runlist))
  fprintf('ERROR: cannot find any runs in %s\n',fsdpath);
  flacnew = [];
  return; 
end

flacnew.nruns = size(flacnew.runlist,1);
if(flacnew.nruns < flac.nthrun)
  fprintf(['ERROR: requested nth run %d is greater than the number of' ...
	   ' runs %d\n'],flacnew.nthrun,flacnew.nruns);
  flacnew = [];
  return; 
end
runid = flacnew.runlist(flac.nthrun,:);

runpath = sprintf('%s/%s',fsdpath,runid);
fstem = sprintf('%s/%s',runpath,flac.funcstem);

% Get the number of time points
mri = MRIread(fstem,1);
if(isempty(mri))
  fprintf('ERROR: attempting to read %s\n',fstem);
  flacnew = [];
  return; 
end
flacnew.mri = mri;
flacnew.ntp = mri.nframes;
flacnew.funcfspec = fstem;

% Time-point exclude file
if(~isempty(flac.tpexcfile))
  fname = sprintf('%s/%s',runpath,flac.tpexcfile);
  fp = fopen(fname,'r'); % ok if it does not exist
  if(fp ~= -1)
    flacnew.tpexc = round(fscanf(fp,'%lf')); 
    flacnew.tpexc = flacnew.tpexc + 1; % change to 1-based
    indtmp = find(flacnew.tpexc >= flacnew.ntp);
    if(~isempty(indtmp))
      fprintf('ERROR: time points in %s exceed nframes (%d)\n',...
	      fname,flacnew.ntp);
      flacnew = [];
      return; 
    end
    fclose(fp);
    % Set up 1-based points to include
    flacnew.tpinc = ones(flacnew.ntp,1); 
    flacnew.tpinc(flacnew.tpexc+1) = 0;
    ntpexc = length(flacnew.tpexc);
    fprintf('Excluding %d points\n',ntpexc);
  else
    flacnew.tpexc = [];
  end
end

% Parfile
if(~isempty(flac.parfile))
  if(isempty(flac.schdir)) % Local
    parpath = sprintf('%s/%s',runpath,flac.parfile);
  else % Global
    parpath = sprintf('%s/r%03d/%s',flac.schdir,flac.nthrun,flac.parfile);
  end
  par = fast_ldpar4(parpath);
  if(isempty(par))
    fprintf('ERROR: loading %s \n',flac.parfile);
    flacnew = []; return; 
  end
  % This is an error check
  if(size(par,2) == 2 & flac.autostimdur)
    % Old-style Two-column par file with autostimdur, check for 0s
    indz = find(par(:,2) == 0);
    if(isempty(indz))
      fprintf('WARNING: %s is a two-column par file. You have specified\n');
      fprintf('that the duration of each stimulus of each presentation\n');
      fprintf('be determined from the par file itself. However, \n');
      fprintf('there are no 0s in this file, so this might not work. \n');
      fprintf('To suppress this msg, consider using a 3 or 4-col par.\n');
    end
  end
end

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  % HRF Regressors
  if(ev.ishrf)  
    % Just get the stimulus timing
    % It is possible that this is empty if the condition is not
    % reprsented for this run. This will cause a bailout unless
    % flac.AllowMissingCond is set.
    if(isempty(flac.parfile))
      % No par file, must be stim timing file
      if(isempty(flac.schdir)) % Local
	stfpath = sprintf('%s/%s',runpath,ev.stf);
      else % Global
	stfpath = sprintf('%s/r%03d/%s',flac.schdir,flac.nthrun,ev.stf);
      end
      st = fast_ldstf(stfpath);
      if(isempty(st) & ~flac.AllowMissingCond)
	fprintf('ERROR: reading timing file %s\n',stfpath);
	flacnew = []; return; 
      end
    else
      % Parfile has been specified
      condno = sscanf(ev.stf,'%d'); % stf is interpreted as cond no
      if(isempty(condno))
	fprintf('ERROR: condition number %s wrong\n',ev.stf);
	flacnew = []; return; 
      end
      if(size(par,2) == 2)
	% Two-column par file
	trun = flacnew.ntp * flacnew.TR;
	if(flac.autostimdur) st = fast_par2st(par,condno,trun);
	else                 st = fast_par2st(par,condno,trun,flac.TR);
	end
      else
	% Three or Four-column par file
	ind = find(par(:,2) == condno);
	st = par(ind,[1 3:end]);
      end
      if(isempty(st) & ~flac.AllowMissingCond)
	fprintf('\nERROR: converting par to st\n');
	fprintf('Could not find condition %d in %s\n\n',condno,parpath);
	flacnew = []; return; 
      end
    end
    if(~isempty(st))
      st(:,1) = st(:,1) + flacnew.stimulusdelay;
      if(strcmp(flac.ev(nthev).model,'fir') & size(st,2) > 2)
	fprintf('Ignoring duration for FIR model\n');
	st(:,3) = zeros(size(st,1),1);
      end
    end
    flacnew.ev(nthev).st = st;
    continue;
  end  
  
  % Non-parametric regressors - load and demean matrix
  if(strcmp(ev.model,'nonpar'))
    nonparpath = sprintf('%s/%s',runpath,ev.nonparname);
    [fspec fstem fmt] = MRIfspec(nonparpath);
    if(~isempty(fspec))
      npmri = MRIread(nonparpath);
      extreg = npmri.vol;
    else
      %fprintf(' ... nonpar trying as a bhdr ...\n');
      extreg = fmri_ldbvolume(nonparpath);
      if(isempty(extreg))
	fprintf('ERROR: loading nonpar reg %s\n',nonparpath);
	flacnew = [];
	return;
      end
    end
    X = fast_vol2mat(extreg);
    if(size(X,1) ~= flacnew.ntp)
      fprintf('ERROR: nonpar time point mismatch %s\n',nonparpath);
      flacnew = [];
      return;
    end
    if(size(X,2) < ev.params(1))
      fprintf('ERROR: not enough columns %s\n',nonparpath);
      size(X)
      keyboard
      flacnew = [];
      return;
    end
    % Demean
    X = X(:,1:ev.params(1));
    Xmn = mean(X,1);
    X = X - repmat(Xmn,[flacnew.ntp 1]);
    flacnew.ev(nthev).X = X;
    continue;
  end

  % Self-Regressor based on Seg
  if(strcmp(ev.model,'selfregseg'))
    if(flac.perrun) runtmp = runid;
    else            runtmp = '';
    end
    segstem = sprintf('%s/%s/masks/%s',fsdpath, runtmp, ev.segstem);
    seg = MRIread(segstem);
    if(isempty(seg))
      fprintf('ERROR: flac_customize: reading %s\n',segstem);
      flacnew = [];
      return;
    end
    segid = ev.params(2:end);
    nseg = length(segid);
    segmask = zeros(size(seg.vol));
    for nthseg = 1:nseg
      segmask = (segmask | seg.vol == segid(nthseg));
    end
    indseg = find(segmask);
    nindseg = length(indseg);
    if(nindseg == 0)
      fprintf('ERROR: flac_customize: %s: no voxels found in seg \n',ev.name);
      flacnew = [];
      return;
    end
    fmri = MRIread(fstem);
    f = fast_vol2mat(fmri.vol);
    clear fmri;
    fseg = f(:,indseg);
    fseg = fseg - repmat(mean(fseg,2),[1 nindseg]);
    [u s v] = fast_svd(fseg);
    npca = ev.params(2);
    if(size(u,2) < npca)
      fprintf('ERROR: flac_customize: %s: not enough components \n',ev.name);
      flacnew = [];
      return;
    end
    X = u(:,1:npca);
    flacnew.ev(nthev).X = X;
    continue;
  end

end

% Now create the full design matrix. This will also create the
% matrices for the HRF-related EVs.
flacnew = flac_desmat(flacnew);
flacnew.indtask = flac_taskregind(flacnew);			    
flacnew.indnuis = flac_nuisregind(flacnew);			    

flacnew.betafspec = sprintf('%s/%s/%s/%s/beta',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));
flacnew.rvarfspec = sprintf('%s/%s/%s/%s/rvar',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));

flacnew.resfspec = sprintf('%s/%s/%s/%s/res',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));

if(~isempty(flacnew.mask))
  flacnew.maskfspec = sprintf('%s/%s/masks/%s%s',flacnew.sess,...
			      flacnew.fsd,flacnew.mask);
else
  flacnew.maskfspec = '';
end
  
flacnew.acfsegfspec = sprintf('%s/%s/masks/%s',flacnew.sess,...
		      flacnew.fsd,flacnew.acfsegstem);

%fprintf('\n');
%fprintf('run %d ---------------------- \n',flacnew.nthrun);
ncon = length(flacnew.con);
for nthcon = 1:ncon
  flacnew.con(nthcon).ffspec = ...
      sprintf('%s/%s/%s/%s/%s/f',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).fsigfspec = ...
      sprintf('%s/%s/%s/%s/%s/fsig',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).gamfspec = ...
      sprintf('%s/%s/%s/%s/%s/ces',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).gamcvmfspec = ...
      sprintf('%s/%s/%s/%s/%s/cesvar',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  % Dont print out efficiency here
  % C = flacnew.con(nthcon).C;
  % eff = 1/trace(C*inv(flacnew.X'*flacnew.X)*C');
  % vrf = 1/mean(diag(C*inv(flacnew.X'*flacnew.X)*C'));
  % flacnew.con(nthcon).eff = eff;
  % fprintf('%2d %s eff=%g, vrf=%g\n',nthcon,flacnew.con(nthcon).name,eff,vrf);
end


return;














