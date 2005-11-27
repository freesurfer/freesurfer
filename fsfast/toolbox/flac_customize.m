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
% $Id: flac_customize.m,v 1.8 2005/11/27 01:44:12 greve Exp $

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

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  % HRF Regressors
  if(ev.ishrf)  
    % Just get the stimulus timing
    stfpath = sprintf('%s/%s',runpath,ev.stf);
    st = fast_ldstf(stfpath);
    if(isempty(st)) 
      fprintf('ERROR: reading timing file %s\n',stfpath);
      flacnew = []; return; 
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

end


flacnew = flac_desmat(flacnew);
flacnew.indtask = flac_taskregind(flacnew);			    
flacnew.indnuis = flac_nuisregind(flacnew);			    

flacnew.betafspec = sprintf('%s/%s/fla/%s/%s/beta',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));
flacnew.rvarfspec = sprintf('%s/%s/fla/%s/%s/rvar',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));

flacnew.resfspec = sprintf('%s/%s/fla/%s/%s/res',flacnew.sess,...
			    flacnew.fsd,flacnew.name,...
			    flacnew.runlist(flacnew.nthrun,:));

flacnew.maskfspec = sprintf('%s/%s/masks/%s%s',flacnew.sess,...
			    flacnew.fsd,flacnew.mask);

flacnew.acfsegfspec = sprintf('%s/%s/masks/%s',flacnew.sess,...
		      flacnew.fsd,flacnew.acfsegstem);

ncon = length(flacnew.con);
for nthcon = 1:ncon
  flacnew.con(nthcon).ffspec = ...
      sprintf('%s/%s/fla/%s/%s/%s/f',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).fsigfspec = ...
      sprintf('%s/%s/fla/%s/%s/%s/fsig',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  flacnew.con(nthcon).gamfspec = ...
      sprintf('%s/%s/fla/%s/%s/%s/gam',flacnew.sess,flacnew.fsd,flacnew.name,...
	      flacnew.runlist(flacnew.nthrun,:),flacnew.con(nthcon).name);
  
end


return;














