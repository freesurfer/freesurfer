% flac_resynth.m
% $Id: flac_resynth.m,v 1.1 2004/12/08 23:36:51 greve Exp $

monly = 1;
flacfile = '~/fbirn-hp-fsfast/flac/sm.flac';
sess = '~/fbirn-hp-fsfast/mgh-data/mgh-103.1';
okfile = '/tmp/flac_resynth.ok';
conmaskname = 'omnibus';
outfspec = 'fmcsm5-sm-resynth';
outconmaskfspec = 'sm-resynth-mask';
fsigthresh = 2;

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

  betafspec = sprintf('%s/%s/fla/%s/%s/beta',...
		      flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  beta = MRIread(betafspec);
  if(isempty(beta)) 
    if(~monly) quit; end
    return;
  end
  beta.vol = fast_vol2mat(beta.vol);

  rvarfspec = sprintf('%s/%s/fla/%s/%s/rvar',...
		      flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  rvar = MRIread(rvarfspec);
  if(isempty(rvar)) 
    if(~monly) quit; end
    return;
  end

  matfile = sprintf('%s/%s/fla/%s/%s/flac.mat',...
		    flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  flacproc = load(matfile);

  
  noise = beta;
  noise.vol = randn(flac.ntp,nv);
  noise.vol = noise.vol .* repmat(sqrt(rvar.vol(:))',[flac.ntp 1]);
  
  nacfseg = size(flacproc.nacfseg,2);
  for nthacfseg = 1:nacfseg
    indseg = find(acfseg.vol==nthacfseg);
    nacf = flacproc.nacfseg(:,nthacfseg);
    S = toeplitz(nacf);
    F = chol(S); % Check this!
    noise.vol(:,indseg) = F*noise.vol(:,indseg);
  end

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
  
end

if(~isempty(outconmaskfspec))
  fname = sprintf('%s/%s/masks/%s%s',flac.sess,flac.fsd,...
		  outconmaskfspec,flac.formatext);
  MRIwrite(conmask,fname);
end

fmri_touch(okfile);

fprintf('flac_resynth: done\n');
