% flac_resynth.m


%
% flac_resynth.m
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

flac_resynth_version = 'flac_resynth.m @FS_VERSION@';

flac_resynth_version_no = sscanf(flac_resynth_version,'%*s %*s %s',1);

%flacfile = '~/fbirn-hp-fsfast/flac/sm.flac';
%sess = '~/fbirn-hp-fsfast/mgh-data/mgh-103.1';
%outfspec = 'fmcsm5-sm-rsyn2';
%outconmaskfspec = 'sm-rsyn-mask';

%flacfile = '~/links/amn/flac/fn.flac';
%sess = '~/links/amn/AMN_01';
%outfspec = 'fmcsm5-fn-rsyn';
%outconmaskfspec = 'fn-rsyn-mask';

%flacfile = '~/links/sg1/workmem/flac/edp.flac';
%sess = '~/links/sg1/workmem/tl20000621';
%outfspec = 'fmcsm5-edp-rsyn';
%outconmaskfspec = 'fn-edp-mask';

%flacfile = '~/links/sg1/xval/flac/samc.flac';
%sess = '~/links/sg1/xval/dng';
%outfspec = 'samc-rsyn';
%outconmaskfspec = 'samc-rsyn-mask';
%noiserlf = 'rest.rlf';

%flacfile = '~/links/sg1/swftst/flac/sm-f00.flac';
%sess = '~/links/sg1/swftst/mgh-103.1';
%outfspec = 'fmc-rsyn';
%outconmaskfspec = 'rsyn-mask';
%noiserlf = 'rest.rlf';

flacfile = '~/links/sg1/mind-swf/flac/sirp.flac';
sess = '~/links/sg1/mind-swf/M02101222';
outfspec = 'fmc-rsyn';
outconmaskfspec = 'rsyn-mask';
noiserlf = 'rest.rlf';

conmaskname = 'probe';
conmap = 'fsigclust.p1';
%conmap = 'fsig';
fsigthresh = 0;
signalscale = 1; % globally scale task betas
noisetempcor = 1; % Temporally correlated noise
noisespatcor = 1; % Spatially  correlated noise
SynthSeed = -1;
nevreg = 10;
usebetaffx = 1;
useconffx = 1;
rvarnorm = 1; % make rest run var same as task run, vox-by-vox
monly = 1;
okfile = '/tmp/flac_resynth.ok';

if(SynthSeed < 0) SynthSeed = round(sum(100*clock)); end
randn('state',SynthSeed); 
rand('state',SynthSeed); 

fprintf('version      %s\n',flac_resynth_version_no);
fprintf('flacfile     %s\n',flacfile);
fprintf('conmask      %s\n',conmaskname);
fprintf('fsigthresh   %g\n',fsigthresh);
fprintf('noiserlf     %s\n',noiserlf);
fprintf('noisetempcor %d\n',noisetempcor);
fprintf('noisespatcor %d\n',noisespatcor);
fprintf('nevreg       %d\n',nevreg);
fprintf('usebetaffx   %d\n',usebetaffx);
fprintf('useconffx    %d\n',useconffx);
fprintf('SynthSeed    %d\n',SynthSeed);

flac = fast_ldflac(flacfile);
if(isempty(flac))
  if(~monly) quit; end
  return;
end

flac.sess = sess;
flac.nthrun = 1;
flac = flac_customize(flac);
if(isempty(flac))
  if(~monly) quit; end
  return;
end
nruns = size(flac.runlist,1);

fsdpath = sprintf('%s/%s',flac.sess,flac.fsd);

% Load rest runs as noise
if(~isempty(noiserlf))
  noiserunlist = fast_runlist(fsdpath,noiserlf);
  if(isempty(noiserunlist))
    fprintf('ERROR: loading %s\n',noiserlf);
    if(~monly) quit; end
    return;
  end
  nnoiseruns = size(noiserunlist,1);
  if(nnoiseruns ~= nruns)
    fprintf('ERROR: number of noise runs (%d) != nruns (%d)\n',...
	    nnoiseruns,nruns);
    if(~monly) quit; end
    return;
  end
  noiseext = 1;
else
  noiseext = 0;
end

%maskfspec = sprintf('%s/bold/masks/synth.mgh',flac.sess);
%mask = MRIread(maskfspec);
mask = MRIread(flac.maskfspec);
if(isempty(mask)) 
  if(~monly) quit; end
  return;
end
indmask = find(mask.vol);

acfseg = MRIread(flac.acfsegfspec);
if(isempty(acfseg)) 
  if(~monly) quit; end
  return;
end

fprintf('Creating contrast mask\n');
conmask = mask;
if(~isempty(conmaskname))
  conindex = flac_conindex(conmaskname,flac);
  if(isempty(conindex))
    if(~monly) quit; end
    return;
  end

  if(~useconffx)
    % Make mask the intersection of each run
    conmask.vol = zeros(size(conmask.vol));
    for nthrun = 1:nruns
      confspec = sprintf('%s/%s/%s/%s/%s/%s',flac.sess,flac.fsd,...
			 flac.name,flac.runlist(nthrun,:),...
			 flac.con(conindex).name,conmap);
      con = MRIread(confspec);
      if(isempty(con)) 
	if(~monly) quit; end
	return;
      end
      conmask.vol = conmask.vol | (abs(con.vol) > fsigthresh);
    end
  else
    % Make mask from FFX
    confspec = sprintf('%s/%s/%s/ffx/%s/%s',flac.sess,flac.fsd,...
		       flac.name,flac.con(conindex).name,conmap);
    con = MRIread(confspec);
    if(isempty(con)) 
      if(~monly) quit; end
      return;
    end
    conmask.vol = abs(con.vol) > fsigthresh;
  end
  
end

conmask.vol = conmask.vol & mask.vol;
nv = prod(conmask.volsize);

indconmask = find(conmask.vol);
indnotconmask = find(~conmask.vol);
nindconmask = length(indconmask);
fprintf('There are %d (%4.1f%%) voxels in the contrast mask\n',...
	nindconmask,100*nindconmask/nv);
if(nindconmask == 0)
  fprintf('ERROR: no point in mask\n');
  if(~monly) quit; end
  return;
end

if(rvarnorm)
  % Load the fixed-effects variance
  rvarfspec = sprintf('%s/%s/ffx/rvar.mgh',fsdpath,...
			 flac.name);
  rvar = MRIread(rvarfspec);
  if(isempty(rvar)) 
    if(~monly) quit; end
    return;
  end
end

%----------------------------------------------------------------%
fprintf('Synthesizing\n');
for nthrun = 1:nruns

  fprintf('nthrun = %d/%d\n',nthrun,nruns);
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  if(isempty(flac))
    if(~monly) quit; end
    return;
  end

  if(noiseext)
    noisefspec = sprintf('%s/%s/%s',fsdpath,...
			 noiserunlist(nthrun,:),flac.funcstem);
    noise = MRIread(noisefspec);
    if(rvarnorm)
      nvar = std(noise.vol,[],4).^2;
      ind = find(nvar==0);
      nvar(ind) = 1000;
      fvar = sqrt(rvar.vol./nvar);
      noise.vol = noise.vol.*repmat(fvar,[1 1 1 flac.ntp]);
    end
    noise.vol = fast_vol2mat(noise.vol);
  else
  
    % Create the noise %
    noise = mask; % Inherit the volume attributes
    noise.vol = randn(flac.ntp,nv);

    % Spatially correlated noise
    if(noisespatcor)
      % Load the resdiual
      resfspec = sprintf('%s/%s/%s/%s/res',...
			 flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
      res = MRIread(resfspec);
      if(isempty(res)) 
	if(~monly) quit; end; 
	return;  
      end
      res.vol = fast_vol2mat(res.vol);
      resmask = res.vol(:,indmask);
      
      % Determine the spatial cov matrix with SVD
      fprintf('Computing SVD of residual\n');
      [Un Sn Vn] = fast_svd(resmask);
      
      % Regularize - saturate at 10th EV
      dSn = diag(Sn);
      dSn(nevreg:end) = dSn(nevreg);
      Sn = diag(dSn);
      %ind = find(dSn < .1*dSn(1));
      %dSn(ind) = .1*dSn(1);
      
      % Spatially filter the white noise
      noise.vol(:,indmask) = ((noise.vol(:,indmask)*Vn)*Sn)*Vn'; % Not Sn.^2
      
    end
    
    
    % Temporally correlated noise
    if(noisetempcor)
      
      matfile = sprintf('%s/%s/%s/%s/flac.mat',...
			flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
      flacproc = load(matfile);
      
      nacfseg = size(flacproc.nacfseg,2);
      for nthacfseg = 1:nacfseg
	indseg = find(acfseg.vol==nthacfseg);
	nacf = flacproc.nacfseg(:,nthacfseg);
	S = toeplitz(nacf);
	F = chol(S); % Check this!
	noise.vol(:,indseg) = F*noise.vol(:,indseg);
      end
    end
    
    % Rescale so that the variance is the same as the original
    rvarfspec = sprintf('%s/%s/%s/%s/rvar',...
			flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
    rvar = MRIread(rvarfspec);
    if(isempty(rvar)) 
      if(~monly) quit; end
      return;
    end
    
    noisestd  = std(noise.vol);
    ind = find(noisestd==0);
    noisestd(ind) = mean(noisestd);
    
    vscale = sqrt(rvar.vol(:))'./noisestd;
    noise.vol = noise.vol .* repmat(vscale,[flac.ntp 1]);
  end % Synthesize Noise

  % ----------- Create the signal -------------------%
  if(usebetaffx)
    betafspec = sprintf('%s/%s/%s/ffx/beta',...
			flac.sess,flac.fsd,flac.name);
  else
    betafspec = sprintf('%s/%s/%s/%s/beta',...
			flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  end
  beta = MRIread(betafspec);
  if(isempty(beta)) 
    if(~monly) quit; end
    return;
  end
  beta.vol = fast_vol2mat(beta.vol);

  % Global scaling of task betas
  beta.vol(flac.indtask,:) = signalscale * beta.vol(flac.indtask,:);
  
  % Comute the signal 
  if(noiseext)
    % Comute the signal in the signal region
    % If using an external source of noise, dont include the
    % nuissance regressors -- it already has nuissance
    sconmask = flac.X(:,flac.indtask)*beta.vol(flac.indtask,indconmask);
    % Comute the signal in the non-signal region - not there, just
    % set to 0
    snotconmask = 0;
  else
    % Purely synthsized noise
    % Comute the signal in the signal region
    sconmask = flac.X*beta.vol(:,indconmask);
    % Comute the signal in the non-signal region as nuissance
    snotconmask = flac.X(:,flac.indnuis)*beta.vol(flac.indnuis,indnotconmask);
  end
  
  y = noise;
  y.vol(:,indconmask)    = y.vol(:,indconmask)    + sconmask; % signal
  y.vol(:,indnotconmask) = y.vol(:,indnotconmask) + snotconmask; % nonsignal

  % Convert back into a volume
  y.vol = fast_mat2vol(y.vol,y.volsize);
  
  outfpath = sprintf('%s/%s/%s/%s%s',flac.sess,flac.fsd,...
		     flac.runlist(nthrun,:),outfspec,flac.formatext);
  MRIwrite(y,outfpath);

  % Write out a logfile too
  outlog = sprintf('%s/%s/%s/%s.log',flac.sess,flac.fsd,...
		   flac.runlist(nthrun,:),outfspec);
  fp = fopen(outlog,'w');
  fprintf(fp,'source       flac_resynth.m\n');
  fprintf(fp,'version      %s\n',flac_resynth_version_no);
  fprintf(fp,'flacfile     %s\n',flacfile);
  fprintf(fp,'conmask      %s\n',conmaskname);
  fprintf(fp,'signalscale  %g\n',signalscale);
  fprintf(fp,'fsigthresh   %g\n',fsigthresh);
  fprintf(fp,'noiserlf     %s\n',noiserlf);
  fprintf(fp,'noisetempcor %d\n',noisetempcor);
  fprintf(fp,'noisespatcor %d\n',noisespatcor);
  fprintf(fp,'nevreg       %d\n',nevreg);
  fprintf(fp,'usebetaffx   %d\n',usebetaffx);
  fprintf(fp,'useconffx    %d\n',useconffx);
  fprintf(fp,'SynthSeed    %d\n',SynthSeed);
  fclose(fp);
  
end % Loop over runs

if(~isempty(outconmaskfspec))
  fname = sprintf('%s/%s/masks/%s%s',flac.sess,flac.fsd,...
		  outconmaskfspec,flac.formatext);
  MRIwrite(conmask,fname);
  fname = sprintf('%s/%s/masks/%s-dil%s',flac.sess,flac.fsd,...
		  outconmaskfspec,flac.formatext);
  conmaskdil1 = conmask;
  conmaskdil1.vol = fast_dilate(conmask.vol,1);
  MRIwrite(conmaskdil1,fname);
end

fmri_touch(okfile);

fprintf('flac_resynth: done\n');
