% fast_stxgrinder2_sess
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
% fast_stxgrinder2_sess.m
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

fprintf('\n\n');
if(~exist('UseMRIread')) UseMRIread=0; end

% Get the output extension
ext = getenv('FSF_OUTPUT_FORMAT');
if(~isempty(ext)) 
  UseMRIread = 1;
else
  ext = 'bhdr'; 
end
fprintf('UseMRIread = %d, ext = %s\n',UseMRIread,ext);


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
    hstem = sprintf('%s/h%s',sessanadir,hemicode);
    h0stem = sprintf('%s/h%s-offset',sessanadir,hemicode);


    % get the dim
    if(~UseMRIread)
      [nrows ncols nframes fs nslices endian bext] = fmri_bfiledim(hstem);
      if(isempty(nrows))
	fprintf('ERROR: loading %s\n',hstem);
	return;
      end
      mristruct = fast_ldbhdr(hstem);
    else
      hmri = MRIread(hstem,1);
      nrows = hmri.volsize(1);
      ncols = hmri.volsize(2);
      nslices = hmri.volsize(3);
      nframes = hmri.nframes;
      mristruct = hmri.bhdr;
    end

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
      J = size(C,1);

      % Compute variance reduction factor
      datfile = sprintf('%s.dat',hstem);
      hd = fmri_lddat3(datfile);
      if(isempty(hd)) return; end
      Ch = hd.hCovMtx;
      concvm = C*Ch*C';
      convrf = 1/mean(diag(concvm));
      fprintf('     VRF: %g\n',convrf);

      if(UseMRIread)
	cesmri = hmri;
	cespctmri  = hmri;
	cesvarmri  = hmri;
	tmri       = hmri;
	sigmri     = hmri;
	minsigmri  = hmri;
	iminsigmri = hmri;
	Fmri       = hmri;
	Fsigmri    = hmri;
      end
      
      
      % ------ Loop over each slice separately ------- %
      fprintf('     slice ');
      for slice = 0:nslices-1
        fprintf('%d ',slice);
	if(rem(slice,21)==20) fprintf('\n           '); end

        % Load beta %
	if(~UseBetaVol)
	  if(~UseMRIread)
	    hAvgFile = sprintf('%s_%03d.bfloat',hstem,slice);
	    [beta rvar hd] = fast_ldsxabfile(hAvgFile);
	    if(isempty(beta))
	      fprintf('ERROR: loading %s\n',hAvgFile);
	      return;
	    end
	  else
	    if(slice == 0) 
	      [hvol rvarvol hd] = fast_ldsxavol(hstem); 
	      if(isempty(hvol)) return; end
	    end
	    %beta = squeeze(hvol.vol(:,:,slice+1,:));
	    % rvar = squeeze(rvarvol.vol(:,:,slice+1));
	    % have to use permute in case there are other singletons
	    beta = permute(hvol.vol(:,:,slice+1,:),[1 2 4 3]); 
	    rvar = permute(rvarvol.vol(:,:,slice+1,:),[1 2 4 3]); 
	  end
	  Ch = hd.hCovMtx;
	  DOF = hd.DOF;
	else
	  betastem = sprintf('%s/beta%s',sessanadir,hemicode);
	  betavarstem = sprintf('%s/beta-var%s',sessanadir,hemicode);
	  beta = fast_ldbslice(betastem,slice);
	  if(isempty(beta))
	    fprintf('ERROR: loading %s\n',betastem);
	    return;
	  end
	  rvar = fast_ldbslice(betavarstem,slice);
	  xmatfile = sprintf('%s/X.mat',sessanadir);
	  XX = load(xmatfile);
	  X = XX.Xfinal;
	  Ch = inv(X'*X);
	  DOF = size(X,1) - size(X,2);
	end
	  
        [nrows ncols nbeta] = size(beta);
        nv = nrows*ncols;
        beta = reshape(beta,[nv nbeta])';
        rvar = reshape(rvar,[nv 1])';
        
	ind = find(rvar == 0 | isnan(rvar));
        rvar(ind) = 10e10;

        % Load mean offset %
	if(~UseMRIread)
	  h0 = fast_ldbslice(h0stem,slice);
	  if(isempty(h0))
	    fprintf('ERROR: loading %s\n',h0stem);
	    return;
	  end
	else
	  if(slice == 0) 
	    hoffsetvol = MRIread(h0stem);
	    if(isempty(hoffsetvol)) return; end
	  end
	  h0 = squeeze(hoffsetvol.vol(:,:,slice+1,:));
	end
        ind = find(h0==0);
        h0(ind) = 10e10;
        h0 = reshape(h0,[nv 1])';

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

        % Go through each row of C separately %
        % Dont have to use a loop here, just easier
        ces    = zeros(J,nv);
        cesvar = zeros(J,nv);
        t      = zeros(J,nv);
        p      = zeros(J,nv);
        for k = 1:J
          Crow = C(k,:);
          cesrow = Crow*beta;
          cesvarrow = rvar * (Crow * Ch * Crow');
          trow = cesrow./sqrt(cesvarrow);
          prow = sign(trow).*tTest(DOF,abs(trow),tTestDOFMax);
          ces(k,:)    = cesrow;
          cesvar(k,:) = cesvarrow;
          t(k,:)      = trow;
          p(k,:)      = prow;
        end

        if(tTestSave & ~strcmp(contrast,'omnibus') & ...
	   ~strcmp(contrast,'zomnibus'))
          cesstem = sprintf('%s/ces%s',condir,hemicode);
          tmp = reshape(ces', [nrows ncols J]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,cesstem,slice,'',mristruct);
	  else
	    cesmri.vol(:,:,slice+1,:) = tmp;
	  end

          cespctstem = sprintf('%s/cespct%s',condir,hemicode);
          tmp = 100*(ces./repmat(h0,[J 1]));
          tmp = reshape(tmp', [nrows ncols J]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,cespctstem,slice,'',mristruct);
	  else
	    cespctmri.vol(:,:,slice+1,:) = tmp;
	  end

          cesvarstem = sprintf('%s/cesvar%s',condir,hemicode);
          tmp = reshape(cesvar', [nrows ncols J]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,cesvarstem,slice,'',mristruct);
	  else
	    cesvarmri.vol(:,:,slice+1,:) = tmp;
	  end

          tstem = sprintf('%s/t%s',condir,hemicode);
          tmp = reshape(t', [nrows ncols J]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,tstem,slice,'',mristruct);
	  else
	    tmri.vol(:,:,slice+1,:) = tmp;
	  end

          pstem = sprintf('%s/sig%s',condir,hemicode);
	  tmp = p; indz = find(p==0); tmp(indz) = 1;
          tmp = -sign(tmp) .* log10(abs(tmp));
          tmp = reshape(tmp', [nrows ncols J]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,pstem,slice,'',mristruct);
	  else
	    sigmri.vol(:,:,slice+1,:) = tmp;
	  end
        end

        % Handle multiple rows in C %
        if(J > 1)
          % Min sig with bonferroni correction
          [ptmp imin] = min(abs(p));
          ind = sub2ind(size(p),imin,1:nv);
          pmin = J*p(ind);
          
          pminstem = sprintf('%s/minsig%s',condir,hemicode);
	  tmp = pmin; indz = find(pmin==0); tmp(indz) = 1;
          tmp = -sign(tmp) .* log10(abs(tmp));
          tmp = reshape(tmp', [nrows ncols 1]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,pminstem,slice,'',mristruct);
	  else
	    minsigmri.vol(:,:,slice+1,:) = tmp;
	  end

          iminstem = sprintf('%s/iminsig%s',condir,hemicode);
          tmp = reshape(imin', [nrows ncols 1]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,iminstem,slice,'bshort',mristruct);
	  else
	    iminsigmri.vol(:,:,slice+1,:) = tmp;
	  end
        end % Handle multiple rows in C

        % F-test
        if(DoFTest)
          cescvm = inv(C*Ch*C');
          if(J>1) F = (sum(ces .* (cescvm*ces))./rvar)/J;
          else    F = t.^2;
          end
	  ind = find(isnan(F));
	  F(ind) = 0;
          Fsig = FTest(J, DOF, F, FTestDOFMax);

          Fstem = sprintf('%s/f%s',condir,hemicode);
          tmp = reshape(F', [nrows ncols 1]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,Fstem,slice,'',mristruct);
	  else
	    Fmri.vol(:,:,slice+1,:) = tmp;
	  end

          Fsigstem = sprintf('%s/fsig%s',condir,hemicode);
	  tmp = Fsig; indz = find(Fsig==0); tmp(indz) = 1;
          tmp = -log10(abs(tmp)); % dont adjust sign
          tmp = reshape(tmp', [nrows ncols 1]);
	  if(~UseMRIread)
	    fast_svbslice(tmp,Fsigstem,slice,'',mristruct);
	  else
	    Fsigmri.vol(:,:,slice+1,:) = tmp;
	  end
        end % FTest
      end % slice
      fprintf('\n');

      if(UseMRIread)
        if(tTestSave & ~strcmp(contrast,'omnibus') & ...
	   ~strcmp(contrast,'zomnibus'))
	  fname = sprintf('%s/ces%s.%s',condir,hemicode,ext);
	  MRIwrite(cesmri,fname);
	  fname = sprintf('%s/cespct%s.%s',condir,hemicode,ext);
	  MRIwrite(cespctmri,fname);
	  fname = sprintf('%s/cesvar%s.%s',condir,hemicode,ext);
	  MRIwrite(cesvarmri,fname);
	  fname = sprintf('%s/t%s.%s',condir,hemicode,ext);
	  MRIwrite(tmri,fname);
	  fname = sprintf('%s/sig%s.%s',condir,hemicode,ext);
	  MRIwrite(sigmri,fname);
	end
	if(J > 1)
	  fname = sprintf('%s/minsig%s.%s',condir,hemicode,ext);
	  MRIwrite(minsigmri,fname);
	  fname = sprintf('%s/iminsig%s.%s',condir,hemicode,ext);
	  MRIwrite(iminsigmri,fname);
	end
        if(DoFTest)
	  fname = sprintf('%s/f%s.%s',condir,hemicode,ext);
	  MRIwrite(Fmri,fname);
	  fname = sprintf('%s/fsig%s.%s',condir,hemicode,ext);
	  MRIwrite(Fsigmri,fname);
	end
      end
      
    end % loop over contrasts      

  end % hemi
  fprintf('\n');
  fprintf('\n');
end % sess

fprintf('\n');
fprintf('done %g\n',toc);
fprintf('\n');

