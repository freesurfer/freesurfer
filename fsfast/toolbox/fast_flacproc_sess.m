% fast_flacproc_sess
%


%
% fast_flacproc_sess.m
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

% Things to fix:
% 1. FIR - does not work for prestim != 0
% 2. FFX computation of gamma covarmtx not correct
% 3. Saving in sxa format

% flacfile = '$flacfile';
% sess = '$sess';
% Synth    = $Synth;
% SynthAR1 = $SynthAR1;
% npca = $npca;
% npcathresh = $npcathresh;
% svres  = $svres;
% svres0 = $svres0;
% svrvar0 = $svrvar0;
% do_fit = $do_fit;
% do_fla = $do_fla;
% do_rfx = $do_rfx;
% do_ffx = $do_ffx;
% monly = $monly;
tic;

fprintf('-----------------------------------------------------\n');
fprintf('Session: %s \n',sess);
fprintf('%g ',clock);fprintf('\n');

if(~isempty(flacfile))
  fprintf('Loading flac %s\n',flacfile);
  flac = fast_ldflac(flacfile);
else
  fprintf('Loading analysis %s\n',analysis);
  flac = fast_ldanaflac(analysis);
end
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

if(~flac.perrun)
  mstem = sprintf('%s/%s/masks/%s%s',flac.sess,flac.fsd,flac.mask);
  acfsegstem = sprintf('%s/%s/masks/%s',flac.sess,flac.fsd,flac.acfsegstem);
end

nruns = size(flac.runlist,1);
if(~do_fla) nruns=0; end

if(~do_fit) fprintf('Not fitting\n'); end

for nthrun = 1:nruns
  flac.nthrun = nthrun;

  fprintf('---------------------------------------------------\n');
  fprintf('Analyzing nthrun %d/%d (%s) %6.1f ------------\n',...
          nthrun,nruns,flac.runlist(nthrun,:),toc);
  fprintf('FLAC: %s, %s\n',flac.name,flac.sess);
  outdir = sprintf('%s/%s/%s/%s',flac.sess,flac.fsd,...
		   flac.name,flac.runlist(flac.nthrun,:));
  mkdirp(outdir);

  flac = flac_customize(flac);
  if(isempty(flac)) 
    if(~monly) quit;  end
    return; 
  end

  if(flac.perrun)
    runid = flac.runlist(flac.nthrun,:);
    mstem = sprintf('%s/%s/%s/masks/%s%s',flac.sess,...
		    flac.fsd,runid,flac.mask);
    acfsegstem = sprintf('%s/%s/%s/masks/%s',flac.sess,flac.fsd, ...
			 runid,flac.acfsegstem);
  end

  if(flac.perrun | nthrun == 1)
    fprintf('%s\n',mstem);
    mask = MRIread(mstem);
    if(isempty(mask))
      if(~monly) quit; end
      return;
    end
    indmask = find(mask.vol);
    indnotmask = find(~mask.vol);
    
    fprintf('%s\n',acfsegstem);
    acfseg = MRIread(acfsegstem);
    if(isempty(acfseg)) 
      if(~monly) quit; end
      return; 
    end
    acfseglist = unique(acfseg.vol(:));
    nacfsegs = length(acfseglist)-1; % Exclude 0
    acfmask = (acfseg.vol ~= 0);
    indacfmask = find(acfmask);
  end

  % Check condition of design matrix after normalizing columns
  % Should check each contrast
  Xtmp = flac.X./repmat(sqrt(sum(flac.X.^2)),[flac.ntp 1]);
  c = cond(Xtmp'*Xtmp);
  fprintf('Design matrix condition: %g\n',c);
  if(c > 1000000)
    fprintf('ERROR: design matrix ill-conditioned\n');
    return;
  end
  ystem = sprintf('%s/%s/%s/%s',flac.sess,flac.fsd,...
		  flac.runlist(flac.nthrun,:),flac.funcstem);
  %fprintf('ystem = %s\n',ystem);

  if(do_fit) %-------------------------------------------------
    fprintf('Loading functional data (%6.1f)\n',toc);
    y = MRIread(ystem);
    if(isempty(y)) 
      if(~monly) quit;  end
      return; 
    end
    szvol = size(y.vol);
    szvol = szvol(1:end-1);
    Nv = prod(szvol);
    y.vol = fast_vol2mat(y.vol);
    
    if(flac.inorm > 0)
      ygmn = mean(reshape1d(y.vol(:,indmask)));
      fprintf('Inorming, ygmn = %g \n',ygmn);
      y.vol = y.vol * (flac.inorm/ygmn);
    end
  
    if(Synth)
      fprintf('Synthesizing input, AR1 = %g \n',SynthAR1);
      y.vol = randn(size(y.vol));
      if(SynthAR1 ~= 0)
	acfsynth = SynthAR1.^[0:flac.ntp-1]';
	Facfsynth = chol(toeplitz(acfsynth))';
	y.vol = Facfsynth*y.vol;
      end
    end
    
    fprintf('Performing OLS estimation   (%6.1f)\n',toc);
    [beta, rvar, vdof, r] = fast_glmfitw(y.vol,flac.X);
    indrvarz = find(rvar==0);
    rvar(indrvarz) = 10e10;

    if(svres0)
      stem = sprintf('%s/res0%s',outdir,flac.formatext);
      tmp = y;
      tmp.vol = fast_mat2vol(r,szvol);
      MRIwrite(tmp,stem);
    end
    
    if(svrvar0)
      stem = sprintf('%s/rvar0%s',outdir,flac.formatext);
      tmp = y;
      tmp.vol = fast_mat2vol(rvar,szvol);
      MRIwrite(tmp,stem);
    end
    
    if(svar1)
      % Save unwhitened ar1
      stem = sprintf('%s/ar1%s',outdir,flac.formatext);
      ar1 = sum(r(1:end-1,:).*r(2:end,:))./sum(r.*r);
      tmp = y;
      tmp.vol = fast_mat2vol(ar1,szvol);
      MRIwrite(tmp,stem);
    end

    racfmri = mask;
    racfmri.vol = zeros(10,Nv);

    fprintf('Computing NACF   (%6.1f)\n',toc);
    % Residual forming matrix 
    R = eye(flac.ntp) - flac.X*inv(flac.X'*flac.X)*flac.X';
    nn = [1:flac.ntp]';
    clear racfseg nacfseg;
    for nthseg = 1:nacfsegs
      %indseg = find(acfseg.vol(indacfmask)==nthseg);
      indseg = find(acfseg.vol==nthseg);
      nperseg(nthseg) = length(indseg);
      racf = fast_acorr(r(:,indseg),[],[],flac.tpexc);
      racfmri.vol(:,indseg) = racf(1:10,:);
      %racfseg(:,nthseg)  = mean(racf(:,indseg),2);
      racfseg(:,nthseg)  = mean(racf,2);
      racfkjw = fast_yacf_kjw(racfseg(1:2,nthseg),R);
      ar1(nthseg) = racfkjw(2);
      nacfseg(:,nthseg) = ar1(nthseg).^(nn-1);
      fprintf('  nthseg = %2d, %5d, rar1 = %6.3f, nar1 = %6.3f\n',...
	      nthseg,nperseg(nthseg),racfseg(2,nthseg),ar1(nthseg));
    end
    clear racf;
    
    if(flac.whiten)
      fprintf('Performing GLS estimation   (%6.1f)\n',toc);
      [beta rvar vdof r] = fast_glmfitw(y.vol,flac.X,nacfseg,acfseg.vol);
    else
      fprintf('NOT Whitening   (%6.1f)\n',toc);
    end
    
    % Replace the zeros with means. The exact value does not really
    % matter, but this allows images to still look ok
    indrvarz = find(rvar==0);
    indrvarnz = find(rvar~=0);
    rvar(indrvarz) = mean(rvar(indrvarnz));

    fprintf('Saving estimation   (%6.1f)\n',toc);
    flacmat = sprintf('%s/flac.mat',outdir);  
    if(str2num(version('-release')) < 14)
      save(flacmat,'flac','racfseg','ar1','nacfseg','acfseg','Synth','SynthAR1');
    else
      save(flacmat,'flac','racfseg','ar1','nacfseg','acfseg','Synth','SynthAR1','-v6');
    end

    stem = sprintf('%s/beta%s',outdir,flac.formatext);
    tmp = y;
    tmp.vol = fast_mat2vol(beta,szvol);
    fprintf('Saving beta to %s\n',stem);
    MRIwrite(tmp,stem);

    % Save the baseline separately
    blevind = flac_evindex(flac,'Baseline');
    if(~isempty(blevind))
      blregind = flac_evregind(flac,blevind);
      tmp = y;
      tmp.vol = fast_mat2vol(beta(blregind,:),szvol);
      stem = sprintf('%s/baseline%s',outdir,flac.formatext);
      MRIwrite(tmp,stem);
    end
    
    % Save residual variance and stddev
    stem = sprintf('%s/rvar%s',outdir,flac.formatext);
    tmp = y;
    tmp.vol = fast_mat2vol(rvar,szvol);
    MRIwrite(tmp,stem);
    stem = sprintf('%s/rstd%s',outdir,flac.formatext);
    tmp.vol = sqrt(tmp.vol);
    MRIwrite(tmp,stem);

    % Save residual ACF
    stem = sprintf('%s/racf%s',outdir,flac.formatext);
    racfmri.vol = fast_mat2vol(racfmri.vol,racfmri.volsize);
    MRIwrite(racfmri,stem);

    % OK, Save residual too
    if(svres)
      stem = sprintf('%s/res%s',outdir,flac.formatext);
      tmp = y;
      tmp.vol = fast_mat2vol(r,szvol);
      MRIwrite(tmp,stem);
    end
    %------ do_fit -----------------------------------------
  else 
    fprintf('  Not fitting\n');
    y = MRIread(ystem,1);
	
    flacmat = sprintf('%s/flac.mat',outdir);  
    tmp = load(flacmat);
    racfseg = tmp.racfseg;
    ar1 = tmp.ar1;
    acfseg = tmp.acfseg;
    nacfseg = tmp.nacfseg;

    stem = sprintf('%s/beta%s',outdir,flac.formatext);
    tmp = MRIread(stem);
    beta = fast_vol2mat(tmp.vol);

    stem = sprintf('%s/rvar%s',outdir,flac.formatext);
    tmp = MRIread(stem);
    rvar = fast_vol2mat(tmp.vol);
  
    szvol = size(tmp.vol);
  end 
  
  fprintf('Computing contrasts   (%6.1f)\n',toc);
  ncontrasts = length(flac.con);
  for nthcon = 1:ncontrasts
    fprintf('  %d  %s   (%6.1f)\n',nthcon,flac.con(nthcon).name,toc);

    if(flac.con(nthcon).varsm > 0)
      fprintf('    Smoothing var with fwhm=%g mm\n',flac.con(nthcon).varsm);
      cfwhm = flac.con(nthcon).varsm/mri.volres(1);
      rfwhm = flac.con(nthcon).varsm/mri.volres(2);
      sfwhm = flac.con(nthcon).varsm/mri.volres(3);
      rvartmp = fast_mat2vol(rvar,szvol);
      rvartmp = fast_smooth3d(rvartmp,cfwhm,rfwhm,sfwhm);
      rvartmp = fast_vol2mat(rvartmp);
    else
      rvartmp = rvar;
    end

    C = flac.con(nthcon).C;
    J = size(C,1);
    if(flac.whiten)
      [F dof1 dof2 ces cescvm] = fast_fratiow(beta,flac.X,rvartmp,C,...
					      nacfseg,acfseg.vol);
    else
      [F dof1 dof2 ces cescvm] = fast_fratiow(beta,flac.X,rvartmp,C);
    end

    p = FTest(dof1, dof2, F);
    sig = -log10(p);

    % If a t-test, then sign it based on the ces
    if(J==1) sig = sig.*sign(ces(1,:)); end

    condir = sprintf('%s/%s',outdir,flac.con(nthcon).name);
    mkdirp(condir);

    stem = sprintf('%s/f%s',condir,flac.formatext);
    tmp = y;
    tmp.vol = fast_mat2vol(F,szvol);
    MRIwrite(tmp,stem);
    if(J==1)
      % For t-test, save as signed t as well as F
      stem = sprintf('%s/t%s',condir,flac.formatext);
      tmp.vol = fast_mat2vol(sqrt(F).*sign(ces(1,:)),szvol);
      MRIwrite(tmp,stem);
    end
    
    
    stem = sprintf('%s/fsig%s',condir,flac.formatext);
    tmp = y;
    tmp.vol = sig;
    if(size(C,1)==1) tmp.vol = tmp.vol .* sign(ces); end
    tmp.vol = fast_mat2vol(tmp.vol,szvol);
    MRIwrite(tmp,stem);

    stem = sprintf('%s/ces%s',condir,flac.formatext);
    tmp = y;
    tmp.vol = fast_mat2vol(ces,szvol);
    MRIwrite(tmp,stem);

    stem = sprintf('%s/cesvar%s',condir,flac.formatext);
    tmp = y;
    tmp.vol = fast_mat2vol(cescvm,szvol);
    MRIwrite(tmp,stem);

    if(dof1 == 1)
      stem = sprintf('%s/cesstd%s',condir,flac.formatext);
      tmp = y;
      tmp.vol = fast_mat2vol(sqrt(cescvm),szvol);
      MRIwrite(tmp,stem);
    end
  
  end % contrasts

  if(npca)
    fprintf('NSVD   (%6.1f)\n',toc);

    % Find task-related voxels
    indtask = flac_taskregind(flac);
    nregtot = size(flac.X,2);
    C = eye(nregtot);
    C = C(indtask,:);
    [F dof1 dof2 ces] = fast_fratiow(beta,flac.X,rvar,C,nacfseg,acfseg.vol);
    p = FTest(dof1, dof2, F);
    sig = -log10(p);
    indnpca = find(abs(sig) < npcathresh);
    fprintf('NSVD: Number of task-related voxels: %d (%g)\n',...
            length(indnpca),100*length(indnpca)/Nv);

    % Extract non-taskrelated voxels and ortho wrt nuissance
    ynpca = y.vol(:,indnpca);
    indnuis = flac_nuisregind(flac);
    Xnuis = flac.X(:,indnuis);
    ynpca = (eye(flac.ntp) - Xnuis*inv(Xnuis'*Xnuis)*Xnuis')*ynpca;
    [u s v] = fast_svd(ynpca);
    ds = diag(s);
    cpvs = 100*cumsum(ds)/sum(ds);
    fprintf('CPVS: ');
    fprintf(' %4.1f',cpvs(1:10));
    fprintf('\n');

    % Keep the all components
    Xnpca = u;
    Xnpca = permute(Xnpca,[3 2 4  1]);
    Xnpcastem = sprintf('%s/Xnpca',outdir);
    fast_svbslice(Xnpca,Xnpcastem);

  end

end % runs

if(do_ffx)
  fprintf('\n\n---------------Computing FFX---------------\n');
  flacffx(flac,'',1);
end

if(do_rfx & nruns > 1)
  fprintf('\n\n---------------Computing RFX---------------\n');
  flacrfx(flac);
end

