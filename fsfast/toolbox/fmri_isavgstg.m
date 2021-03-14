% fmri_isavgstg (Inter-Subject Averaging, Random Effects)
%
%


%
% fmri_isavgstg.m
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

fprintf('\n$Id: fmri_isavgstg.m,v 1.4 2011/03/02 00:04:06 nicks Exp $\n');

nS = size(isavgdat,1);
tmp = load(CMtxFile);
RM = tmp.ContrastMtx_0;

% ---------- Go through each slice ----------------- %
for slc = FirstSlice : FirstSlice+nSlices-1,
  fprintf(1,'slc = %2d  --------------- \n',slc);
  vall = [];

  % ----------- Load slice data for each subject ---------- %
  for s = 1:nS
    fprintf(1,'   session = %2d \n',s);

    if(strcmp(VolType,'selxavg') | strcmp(VolType,'selavg'))
      % ------- Load the .dat file ---------- %
      datfile = sprintf('%s.dat',isavgdat(1).hdrstem);
      hdrdat = fmri_lddat3(datfile);
      Nc = hdrdat.Nc;
      Nh = hdrdat.Nh;
      Nch = (Nc-1)*Nh;
    
      if(Nch ~= size(RM,2))
        msg = sprintf('ERROR: Nch is inconsistent (%d,%d)',Nch,size(RM,2));
        qoe(msg);error(msg);
      end

      % ---- Load data --------%
      hsafile = sprintf('%s_%03d.bfloat',isavgdat(s).hdrstem,slc);
      hsa = fmri_ldbfile(hsafile);

      % ------- Extract averages and stddevs -----------%
      [havg hstd hdrdat] = fmri_untangle(hsa,hdrdat);

      % ------- Subjtract Condition 0 ----------%
      havg0 = havg(:,:,1,:);
      havg0r = repmat(havg0,[1 1 Nc-1 1]);
      havg = havg(:,:,[2:Nc],:) - havg0r;
      Nc = Nc-1;
      clear havg0 havg0r;
    else
      hsafile = sprintf('%s_%03d.bfloat',isavgdat(s).hdrstem,slc);
      havg = fmri_ldbfile(hsafile);
      hdrdat = fmri_hdrdatstruct;
      hdrdat.Nrows = size(havg,1);
      hdrdat.Ncols = size(havg,2);
      Nch = size(havg,3);
    end

    % -------- Reshape ---------------%
    Nv = hdrdat.Nrows * hdrdat.Ncols;
    havg = permute(havg, [4 3 1 2]);
    havg = reshape(havg, [Nch Nv]);

    %------- Synthesize Data ------%
    if(showfpr)
      havg = randn(size(havg));
    end

    if(~isempty(inweights))
      havg = havg*inweights(s);
    end

    % ------ Project using Restriction Matrix -------- %
    v = RM*havg;
    vall = [vall; v];

    if(jackknife)
      if(s==1 & slc == FirstSlice)  
         Rh_all = zeros(size(RM,1),Nv,nS); 
      end
      Rh_all(:,:,s) = fmri_norm(v,2);
    end

  end % end for S
  % ---------- Finished loading data ------------------ %
  
  % ---------- Jackknifing ------------------ %
  if(jackknife)
    [vavg vstd] = fmri_jackknife(vall);
  else
    vavg = mean(vall);
    vstd = std(vall);
    %vstd = sqrt(mean(vall.^2));
  end
  dof = size(vall,1) - 1;

  % ------ check for any data equalling zero -------- %
  ind = find(vstd == 0);
  vstd(ind) = 10.^(-10);

  % -------- compute t values ---------- %
  t = vavg./(vstd/sqrt(dof));
  fprintf('std(t) = %g\n',std(reshape1d(t)));

  % ------- compute p values ----------- %
  p = tTest(dof, t, 100);
  ind = find(p == 0);
  p(ind) = 1;

  if(showfpr)
    h1 = figure(1);
    [FPR alpha] = ComputeFPR(reshape1d(p));
    alphamax = .3;
    nshow = find(alpha<alphamax);
    plot(alpha(nshow),FPR(nshow),'r',alpha(nshow),alpha(nshow),'g+-')
    if(~monly)
      uiwait(h1);
      quit;
    end
    return;
  end

  % ------- save results -------- %
  glbavg = reshape(vavg, [hdrdat.Nrows hdrdat.Ncols]);
  glbstd = reshape(vstd, [hdrdat.Nrows hdrdat.Ncols]);
  t = reshape(t, [hdrdat.Nrows hdrdat.Ncols]);
  p = reshape(p, [hdrdat.Nrows hdrdat.Ncols]);
  if(InvertSig) p = sign(p).*(1-abs(p)); end
  if(strcmp(StatFormat,'ln'))
    fprintf('INFO: converting sig values to natural log\n');
    p = - sign(t) .* log(abs(p));
  elseif(strcmp(StatFormat,'log10'))
    fprintf('INFO: converting sig values to log10\n');
    p = - sign(t) .* log10(abs(p));
  else
    p = - sign(t) .* abs(p);
  end

  OutFile = sprintf('%s_%03d.bfloat',OutStem,slc);
  fmri_svbfile(p,OutFile);

  OutFile = sprintf('%s-avgci_%03d.bfloat',OutStem,slc);
  tmp(:,:,1) = glbavg;
  tmp(:,:,2) = glbstd;
  tmp(:,:,3) = t;
  fmri_svbfile(tmp,OutFile);

end % for slc

fprintf('fmri_isavgstg: done  \n');
if(~monly)
  quit;
end
