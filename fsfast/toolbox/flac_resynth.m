% flac_resynth.m

flac_resynth_version = '$Id: flac_resynth.m,v 1.3 2004/12/12 02:20:15 greve Exp $';

flac_resynth_version_no = sscanf(flac_resynth_version,'%*s %*s %s',1);

%flacfile = '~/fbirn-hp-fsfast/flac/sm.flac';
%sess = '~/fbirn-hp-fsfast/mgh-data/mgh-103.1';
%outfspec = 'fmcsm5-sm-rsyn';
%outconmaskfspec = 'sm-rsyn-mask';

%flacfile = '~/links/amn/flac/fn.flac';
%sess = '~/links/amn/AMN_01';
%outfspec = 'fmcsm5-fn-rsyn';
%outconmaskfspec = 'fn-rsyn-mask';

flacfile = '~/links/sg1/workmem/flac/edp.flac';
sess = '~/links/sg1/workmem/tl20000621';
outfspec = 'fmcsm5-edp-rsyn';
outconmaskfspec = 'fn-edp-mask';

conmaskname = 'omnibus';
fsigthresh = 2;
noisetempcor = 1; % Temporally correlated noise
noisespatcor = 1; % Spatially  correlated noise
SynthSeed = -1;
nevreg = 10;
monly = 1;
okfile = '/tmp/flac_resynth.ok';

if(SynthSeed < 0) SynthSeed = round(sum(100*clock)); end
randn('state',SynthSeed); 
rand('state',SynthSeed); 

fprintf('version      %s\n',flac_resynth_version_no);
fprintf('flacfile     %s\n',flacfile);
fprintf('conmask      %s\n',conmaskname);
fprintf('fsigthresh   %g\n',fsigthresh);
fprintf('noisetempcor %d\n',noisetempcor);
fprintf('noisespatcor %d\n',noisespatcor);
fprintf('nevreg       %d\n',nevreg);
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

flac = flac_desmat(flac);
if(isempty(flac))
  if(~monly) quit; end
  return;
end

maskfspec = sprintf('%s/%s/masks/%s',flac.sess,flac.fsd,flac.mask);
mask = MRIread(maskfspec);
if(isempty(mask)) 
  if(~monly) quit; end
  return;
end
indmask = find(mask.vol);

acfsegfspec = sprintf('%s/%s/masks/%s',flac.sess,flac.fsd,flac.acfsegstem);
acfseg = MRIread(acfsegfspec);
if(isempty(acfseg)) 
  if(~monly) quit; end
  return;
end

fprintf('Creating contrast mask\n');
conmask = mask;
conmask.vol = zeros(size(conmask.vol));
if(~isempty(conmaskname))
  conindex = flac_conindex(conmaskname,flac);
  if(isempty(conindex))
    if(~monly) quit; end
    return;
  end
  for nthrun = 1:nruns
    confspec = sprintf('%s/%s/fla/%s/%s/%s/fsig',flac.sess,flac.fsd,...
		       flac.name,flac.runlist(nthrun,:),...
		       flac.con(conindex).name);
    con = MRIread(confspec);
    if(isempty(con)) 
      if(~monly) quit; end
      return;
    end
    conmask.vol = conmask.vol | (abs(con.vol) > fsigthresh);
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
  

fprintf('Synthesizing\n');
for nthrun = 1:nruns

  fprintf('nthrun = %d/%d\n',nthrun,nruns);
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  if(isempty(flac))
    if(~monly) quit; end
    return;
  end
  flac = flac_desmat(flac);
  if(isempty(flac))
    if(~monly) quit; end
    return;
  end

  % Create the noise %
  noise = mask; % Inherit the volume attributes
  noise.vol = randn(flac.ntp,nv);

  % Spatially correlated noise
  if(noisespatcor)
    % Load the resdiual
    resfspec = sprintf('%s/%s/fla/%s/%s/res',...
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
  
    matfile = sprintf('%s/%s/fla/%s/%s/flac.mat',...
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
  rvarfspec = sprintf('%s/%s/fla/%s/%s/rvar',...
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

  % ----------- Create the signal -------------------%
  betafspec = sprintf('%s/%s/fla/%s/%s/beta',...
		      flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  beta = MRIread(betafspec);
  if(isempty(beta)) 
    if(~monly) quit; end
    return;
  end
  beta.vol = fast_vol2mat(beta.vol);

  sconmask = flac.X*beta.vol(:,indconmask);
  
  indnuis = flac_nuisregind(flac);
  snotconmask = flac.X(:,indnuis)*beta.vol(indnuis,indnotconmask);
  
  y = noise;
  y.vol(:,indconmask) = y.vol(:,indconmask) + sconmask;
  y.vol(:,indnotconmask) = y.vol(:,indnotconmask) + snotconmask;

  y.vol = fast_mat2vol(y.vol,y.volsize);
  
  outfpath = sprintf('%s/%s/%s/%s%s',flac.sess,flac.fsd,...
		     flac.runlist(nthrun,:),outfspec,flac.formatext);
  MRIwrite(y,outfpath);
  % Write out a logfile
  outlog = sprintf('%s/%s/%s/%s.log',flac.sess,flac.fsd,...
		   flac.runlist(nthrun,:),outfspec);
  fp = fopen(outlog,'w');
  fprintf(fp,'source       flac_resynth.m\n');
  fprintf(fp,'version      %s\n',flac_resynth_version_no);
  fprintf(fp,'flacfile     %s\n',flacfile);
  fprintf(fp,'conmask      %s\n',conmaskname);
  fprintf(fp,'fsigthresh   %g\n',fsigthresh);
  fprintf(fp,'noisetempcor %d\n',noisetempcor);
  fprintf(fp,'noisespatcor %d\n',noisespatcor);
  fprintf(fp,'nevreg       %d\n',nevreg);
  fprintf(fp,'SynthSeed    %d\n',SynthSeed);
  fclose(fp);
  
end

if(~isempty(outconmaskfspec))
  fname = sprintf('%s/%s/masks/%s%s',flac.sess,flac.fsd,...
		  outconmaskfspec,flac.formatext);
  MRIwrite(conmask,fname);
end

fmri_touch(okfile);

fprintf('flac_resynth: done\n');
