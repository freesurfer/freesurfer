% fast_fcreg.m
% $Id: fast_fcreg.m,v 1.6 2009/12/14 18:47:32 greve Exp $

%flacfile
%sess 
%segfile
% contrast
% map
% thresh: val, sign
% Left: thal hippo caudmidfront pericalc
% segidlist = [10 17 1003 1021];
% fcregstem = 'fcreg';
% perrun = 1;
% nnuis = 10;


flac = fast_ldflac(flacfile);
if(isempty(flac)) return; end

flac.sess = sess;
flac.nthrun = 1;
flac = flac_customize(flac);
if(isempty(flac)) return; end

nruns = size(flac.runlist,1);
fsdp = sprintf('%s/%s',sess,flac.fsd);

fprintf('\n');
for nthrun = 1:nruns
  runid = flac.runlist(nthrun,:);
  fprintf('%s  run %d %s ---------------- \n',sess,nthrun,runid);
  
  flac.nthrun = nthrun;
  runflac = flac_customize(flac);
  if(isempty(runflac)) return; end
  X = runflac.X;
  nframes = size(X,1);
  R = eye(nframes) - X*inv(X'*X)*X';

  % Load in the segmentation
  if(perrun) 
    segpath = sprintf('%s/%s/masks/%s',fsdp,runid,segfile);
  else
    segpath = sprintf('%s/masks/%s',fsdp,segfile);
  end
  seg = MRIread(segpath);
  if(isempty(seg)) return; end

  fprintf('nSVD = %d\n',nSVD);
  
  if(WMOrthog)
    wm   = (seg.vol ==  2 | seg.vol == 41);
    wm = fast_dilate(wm,1,1,1); % erode 1 in-plane
    indwm = find(wm);
    nwm = length(indwm);
    fprintf('nwm = %d\n',nwm);
  else
    indwm = [];
  end

  if(VCSFOrthog)
    vcsf = (seg.vol ==  4 | seg.vol == 43 | ...
	    seg.vol ==  5 | seg.vol == 44 | ...
	    seg.vol == 14 | seg.vol == 72);
    vcsf = fast_dilate(vcsf,1,1,1); % erode 1 in-plane
    indvcsf = find(vcsf);
    nvcsf = length(indvcsf);
    fprintf('nvcsf = %d\n',nvcsf);
  else
    indvcsf = [];
  end

  f = MRIread(runflac.funcfspec);
  if(isempty(f)) return; end
  fmat  = fast_vol2mat(f);

  Xd = X;
  if(WMOrthog | VCSFOrthog)
    % Construct the nuisance matrix
    rwm   = R*fmat(:,indwm);
    rvcsf = R*fmat(:,indvcsf);
    A = [rwm rvcsf];
    [u s v] = fast_svd(A);
    ds2 = diag(s).^2;
    pvsn = 100*ds2/sum(ds2);
    cpvsn = cumsum(pvsn);
    fprintf('NPVS (%d):  ',nnuis);
    fprintf('%4.1f ',pvsn(1:10));
    fprintf('\n');
    fprintf('NCPVS (%d):  ',nnuis);
    fprintf('%4.1f ',cpvsn(1:10));
    fprintf('\n');

    Xd = [Xd u(:,[1:nnuis])];
  end
  Rd = eye(nframes) - Xd*inv(Xd'*Xd)*Xd';
  
  fcreg = [];
  nthseg = 0;
  for segid = segidlist
    nthseg = nthseg + 1;
    indseg = find(seg.vol == segid);
    nseg = length(indseg);
    if(nseg == 0)
      fprintf('WARNING: cannot find any segid = %d\n',segid);
      %return;
      continue;
    end
    fprintf('%2d segid = %d,  nseg = %g\n',nthseg,segid,nseg);
    fseg = fmat(:,indseg);
    if(nSVD == 0) 
      % Take the mean over the given Seg
      rseg = mean(Rd*fseg,2);
      rseg = rseg/sqrt(sum(rseg.^2));
      fcreg = [fcreg rseg];
    else
      % Concat all waveforms in the Seg for later SVD
      fcreg = [fcreg Rd*fseg];
    end
  end

  if(nSVD ~= 0) 
    [Utmp Stmp Vtmp] = svd(fcreg);
    pvs = 100*diag(Stmp)/sum(diag(Stmp));
    cpvs = cumsum(pvs);
    fprintf('SVD CPVS: %4.1f %4.1f %4.1f %4.1f %4.1f\n',cpvs(1:5));
    fprintf('SVD CPVS Tot (n=%d): %4.1f\n',nSVD,cpvs(nSVD));
    fcreg = Utmp(:,1:nSVD);
  end
  
  fname = sprintf('%s/%s/%s',fsdp,runid,fcregstem);
  fast_svextreg(fcreg,fname);
  
  fprintf('\n');
end
if(exist('okfile','var'))  fmri_touch(okfile); end
















