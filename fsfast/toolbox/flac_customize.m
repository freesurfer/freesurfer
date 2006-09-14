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
% $Id: flac_customize.m,v 1.14 2006/09/14 02:02:55 greve Exp $

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
    

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  % HRF Regressors
  if(ev.ishrf)  
    % Just get the stimulus timing
    if(isempty(flac.parfile))
      stfpath = sprintf('%s/%s',runpath,ev.stf);
      st = fast_ldstf(stfpath);
      if(isempty(st)) 
	fprintf('ERROR: reading timing file %s\n',stfpath);
	flacnew = []; return; 
      end
    else
      parpath = sprintf('%s/%s',runpath,flac.parfile);
      par = fmri_ldpar(parpath);
      if(isempty(par))
	fprintf('ERROR: loading %s \n',flac.parfile);
	flacnew = []; return; 
      end
      condno = sscanf(ev.stf,'%d');
      if(isempty(condno))
	fprintf('ERROR: condition number %s wrong\n',ev.stf);
	flacnew = []; return; 
      end
      trun = flacnew.ntp * flacnew.TR;
      st = fast_par2st(par,condno,trun);
      if(isempty(st))
	fprintf('ERROR: converting par to st\n');
	flacnew = []; return; 
      end
    end
    st(:,1) = st(:,1) + flacnew.stimulusdelay;
    flacnew.ev(nthev).st = st;
    continue;
  end  
  
  % Non-parametric regressors - load and demean matrix
  if(strcmp(ev.model,'nonpar'))
    nonparpath = sprintf('%s/%s',runpath,ev.nonparname);
    X = fast_ldbslice(nonparpath);
    if(isempty(X))
      fprintf('ERROR: loading %s\n',nonparpath);
      flacnew = [];
      return;
    end
    if(size(X,2)~=1) X = squeeze(X)';
    else             X = squeeze(X);
    end
    if(size(X,1) ~= flacnew.ntp)
      fprintf('ERROR: nonpar time point mismatch %s\n',nonparpath);
      flacnew = [];
      return;
    end
    if(size(X,2) < ev.params(1))
      fprintf('ERROR: not enough columns %s\n',nonparpath);
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

flacnew.maskfspec = sprintf('%s/%s/masks/%s%s',flacnew.sess,...
			    flacnew.fsd,flacnew.mask);

flacnew.acfsegfspec = sprintf('%s/%s/masks/%s',flacnew.sess,...
		      flacnew.fsd,flacnew.acfsegstem);

fprintf('\n');
fprintf('run %d ---------------------- \n',flacnew.nthrun);
ncon = length(flacnew.con);
for nthcon = 1:ncon
  flacnew.con(nthcon).ffspec = ...
      sprintf('%s/%s/%s/%s/%s/f',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).fsigfspec = ...
      sprintf('%s/%s/%s/%s/%s/fsig',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).gamfspec = ...
      sprintf('%s/%s/%s/%s/%s/gam',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  C = flacnew.con(nthcon).C;
  eff = 1/trace(C*inv(flacnew.X'*flacnew.X)*C');
  flacnew.con(nthcon).eff = eff;
  fprintf('%2d %s eff=%g\n',nthcon,flacnew.con(nthcon).name,eff);
end


return;














