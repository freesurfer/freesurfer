% fast_stxgrinder3_sess
%
%
% Creates anadir/contrast
%     Univariate:       sig,t,ces,cesvar; sig has same sign as ces
%     Multivariate:     fsig,F,minsig,iminsig
%
% These variables must be defined previously
% SessList = splitstring('$SessList');
% fsd      = '$fsd';
% analysis = '$analysis';
% contrasts = splitstring('$contrastlist');
% hemi = splitstring('$hemi');
% spacedir = '$spacedir';
% tTestDOFMax = $tTestDOFMax;
% FTestDOFMax = $FTestDOFMax;
% DoFTest = $DoFTest;
% tTestSave = $tTestSave;
% IsGroup = [$IsGroupList];
% UseBetaVol = 1;
% OutDir = [];


%
% fast_stxgrinder3_sess.m
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

tic;
nsess = size(SessList,1);
nhemi = size(hemi,1);
ncontrasts = size(contrasts,1);

% Get the output extension
ext = getenv('FSF_OUTPUT_FORMAT');
if(~isempty(ext)) 
  UseMRIread = 1;
else
  ext = 'bhdr'; 
end

for nthsess = 1:nsess
  sessdir = deblank(SessList(nthsess,:));
  sessid = basename(sessdir);
  fprintf('nthsess = %d  %s time=%g --------\n',nthsess,sessid,toc);
  fprintf('%s\n',sessdir);

  for nthhemi = 1:nhemi

    hid = deblank(hemi(nthhemi,:));
    if(strcmp(hid,'nohemi'))  
      hemicode = '';
    else                       
      fprintf('hemi = %s   (%g)\n',hid,toc);
      hemicode = sprintf('-%s',hid);
    end

    if(IsGroup(nthsess))
      sessanadir = sprintf('%s/%s/%s/%s-ffx',sessdir,fsd,analysis,spacedir);
    else
      sessanadir = sprintf('%s/%s/%s/%s',sessdir,fsd,analysis,spacedir);
    end
    xfile = sprintf('%s/X.mat',sessanadir);
    XX = load(xfile);
    if(isempty(XX))
      fprintf('ERROR: could not load %s\n',xfile);
      return;
    end
    X = XX.Xfinal;
    nTask = XX.Navgs_per_cond * XX.Nnnc;
    nNuis = size(X,2) - nTask;
    fprintf('nTask = %d, nNuis = %d\n',nTask,nNuis);

    hstem = sprintf('%s/h%s',sessanadir,hemicode);
    h0stem = sprintf('%s/h%s-offset',sessanadir,hemicode);

    h0 = MRIread(h0stem);
    if(isempty(h0)) return; end
    
    [beta rvar hd] = fast_ldsxavol(hstem); 
    if(isempty(beta)) return; end
    rvar.volmat = fast_vol2mat(rvar.vol);
    beta.volmat = fast_vol2mat(beta.vol);
    [nbeta nvox] = size(beta.volmat);
    beta.volmat = [beta.volmat; zeros(nNuis,nvox)]; % Add 0s for nuis

    acfmatfile = sprintf('%s/acf.mat',sessanadir);
    if(fast_fileexists(acfmatfile))
      fprintf('Loading ACF/AR1\n');
      acfmat = load(acfmatfile);
      if(isempty(acfmat))
	fprintf('ERROR: Loading %s\n',acfmatfile);
	return;
      end
      nacf = length(acfmat.ar1seg);
      nX = size(X,1);
      nn = repmat([1:nX]',[1 nacf]);
      rho = repmat(acfmat.ar1seg,[nX 1]);
      acf = rho.^(nn-1);
      acfseg = acfmat.acfseg.vol;
    else
      fprintf('Not Loading ACF/AR1\n');
      acf = [];
      acfseg = [];
    end
    
    % Initialize
    ces     = beta;
    cesvar  = beta;
    p       = beta;
    sig     = beta;
    minsig  = beta;
    iminsig = beta;
    F       = beta;
    Fsig    = beta;

    % Compute variance reduction factor
    Ch = hd.hCovMtx;

    % Contrast Loop
    for c = 1:ncontrasts
      contrast = deblank(contrasts(c,:));
      fprintf('  contrast %s  (time=%g)\n',contrast,toc);
      if(isempty(OutDir))
	condir = sprintf('%s/%s',sessanadir,contrast);
      else
	condir = sprintf('%s/%s/%s/%s/%s',...
			 OutDir,sessid,fsd,analysis,contrast);
	fprintf('condir %s\n',condir);
      end
      
      cmat = sprintf('%s/%s.mat',analysis,contrast);
      tmp = load(cmat);
      if(isempty(tmp))
        fprintf('ERROR: loading %s\n',cmat);
        return;
      end
      C = tmp.ContrastMtx_0;
      if(size(C,2) ~= nbeta)
	fprintf('\n');
	fprintf('ERROR: size mismatch between analysis %s and contrast %s.\n',...
		analysis,contrast);
	fprintf('This usually happens when the parameters of an analysis\n');
	fprintf('have been changed without re-creating the contrast,\n');
	fprintf('or the analysis was changed and the contrast updated\n');
	fprintf('but selxavg was not re-run for this subject.\n');
	fprintf('\n');
	fprintf('Try re-running mkcontrast-sess for this contrast\n');
	fprintf('and/or re-running selxavg-sess for this subject.\n');
	fprintf('\n');
	return;
      end      
      J = size(C,1);
      C = [C zeros(J,nNuis)]; % Add 0s for nuis

      [F.volmat dof1 dof2 ces.volmat cesvar.volmat] = ...
	  fast_fratiow(beta.volmat,X,rvar.volmat,C,acf,acfseg);
      p.volmat = FTest(dof1, dof2, F.volmat);
      sig.vol = -log10(fast_mat2vol(p));
      if(J > 1)
	% Multivariate test
	fname = sprintf('%s/F%s.%s',condir,hemicode,ext);
	F.vol = fast_mat2vol(F);
	MRIwrite(F,fname);
	fname = sprintf('%s/fsig%s.%s',condir,hemicode,ext);
	MRIwrite(sig,fname);
	p.volmat = [];
	ces.volmat = [];
	for k = 1:J
	  Crow = C(k,:);
	  [F.volmat dof1 dof2 cestmp] = ...
	      fast_fratiow(beta.volmat,X,rvar.volmat,Crow);
	  ptmp = FTest(dof1, dof2, F.volmat);
	  p.volmat = [p.volmat; ptmp];
	  ces.volmat = [ces.volmat; cestmp];
	end
	p.vol = fast_mat2vol(p);
	sig.vol = -log10(p.vol);
	fname = sprintf('%s/sig%s.%s',condir,hemicode,ext);
	MRIwrite(sig,fname);
	
	pmin  = sig;
	ipmin = sig;
	[pmin.volmat ipmin.volmat] = min(p.volmat);
	pmin.volmat = pmin.volmat*J; % Bonf
	indpmin = sub2ind(size(p.volmat),ipmin.volmat,1:nvox);
	cespmin = ces.volmat(indpmin);
	minsig = sig;
	minsig.volmat = -log10(pmin.volmat) .* sign(cespmin);
	minsig.vol = fast_mat2vol(minsig);
	fname = sprintf('%s/minsig%s.%s',condir,hemicode,ext);
	MRIwrite(minsig,fname); % Need to sign minsig
	fname = sprintf('%s/iminsig%s.%s',condir,hemicode,ext);
	ipmin.vol = fast_mat2vol(ipmin);
	MRIwrite(ipmin,fname);
      else
	% T-test
	fname = sprintf('%s/ces%s.%s',condir,hemicode,ext);
	ces.vol = fast_mat2vol(ces);
	MRIwrite(ces,fname);
	fname = sprintf('%s/cespct%s.%s',condir,hemicode,ext);
	cespct = ces;
	cespct.vol = 100*ces.vol./h0.vol;
	MRIwrite(cespct,fname);
	fname = sprintf('%s/cesvar%s.%s',condir,hemicode,ext);
	cesvar.vol = fast_mat2vol(cesvar);
	MRIwrite(cesvar,fname);
	fname = sprintf('%s/cesvarpct%s.%s',condir,hemicode,ext);
	cesvarpct = cesvar;
	cesvarpct.vol = 100*cesvar.vol./h0.vol;
	MRIwrite(cesvarpct,fname);
	t=F;
	t.vol = sqrt(fast_mat2vol(F)) .* sign(ces.vol);
	fname = sprintf('%s/t%s.%s',condir,hemicode,ext);
	MRIwrite(t,fname);
	sig.vol = sig.vol .* sign(ces.vol);
	fname = sprintf('%s/sig%s.%s',condir,hemicode,ext);
	MRIwrite(sig,fname);
      end
    end % loop over contrasts      

  end % hemi
  fprintf('\n');
  fprintf('\n');
end % sess

fprintf('\n');
fprintf('done %g\n',toc);
fprintf('\n');

