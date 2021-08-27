% stxgpredict.m (see stxpredict)
% computes a prediction of the statistical map
% given the averages from an already analyzed data set as well
% as a list of paradigm files.

% These must be previously defined:
% hstem
% parlist
% contrast
% tstem       
% tsigstem    
% tminsigstem 
% fstem       
% fsigstem    
% ntrs - number of time points per run
% firstslice
% nslices


%
% stxgpredict.m
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

% Load the datfile from the analysis %
datfile = sprintf('%s.dat',hstem);
hdr = fmri_lddat3(datfile);

% Number of runs from the list of paradigm files %
nruns = size(parlist,1);

% Make sure that the number of conditions in the 
% paradigm files equals that from the analysis 
parall = [];
for nthrun = 1:nruns
  parfile = deblank(parlist(nthrun,:));
  par = fmri_ldpar(parfile);
  parall = [parall; par];  
end
condlist = unique(par(:,2));
inot0 = find(condlist ~= 0);
condlist = condlist(inot0);
if(length(condlist) ~= hdr.Nnnc)
  fprintf('ERROR: number of conditions in paradigm files (%d)\n',...
          length(condlist));
  fprintf('Does not match that of the analyzed data (%d)\n',hdr.Nnnc);
  return;
end
  
SubSampRate = hdr.TR/hdr.TER;
Nfir = hdr.TimeWindow/hdr.TER;

% Determine whether the mean and trend were fit %
if(hdr.DTOrder == hdr.Nruns) 
   FitMean = 1;
elseif(hdr.DTOrder == 2*hdr.Nruns) 
   FitMean  = 1;
   FitTrend = 1;
else
   FitMean  = 0;
   FitTrend = 0;
end

% Go through each run and compute the Xfir and detrending matrices
Xfir = [];
Xdt = [];
for nthrun = 1:nruns
  parfile = deblank(parlist(nthrun,:));
  par = fmri_ldpar(parfile);
  Xfirrun = fmri_par2scm(par,hdr.Nc,SubSampRate*ntrs,hdr.TER,...
    Nfir,hdr.TPreStim);
  Xfir = [Xfir; Xfirrun];

  Xbaseline = []; Xtrend    = [];
  if(FitMean)  Xbaseline = fast_baselinemtx(nthrun,ntrs,nruns); end
  if(FitTrend) Xtrend    = fast_trendmtx(nthrun,ntrs,nruns); end
  Xdt = [Xdt; [Xbaseline Xtrend]];
end

% For Sub-TR Estimation %
if(hdr.TR ~= hdr.TER)
  Xfirtmp = Xfir;
  nn = [1:SubSampRate:size(Xfirtmp,1)];
  Xfir = Xfirtmp(nn,:);
end

% Final Desgin Matrix %
X = [Xfir Xdt];

% load contrast %
s = load(contrast);
R = s.ContrastMtx_0;

XtX = X'*X; %'
iXtX = inv(XtX);
nkeep = size(iXtX,1) - nruns*(FitMean + FitTrend);
iXtX = iXtX(1:nkeep,1:nkeep);

% Check consistency %
if(size(R,2) ~= nkeep)
  fprintf('ERROR: size(R,2) ~= nkeep (%d,%d)\n',size(R,2),nkeep);
  return;
end

% Compute CF which makes the calucation of the F-stat easier:
% F = (CF*R*h)' * (CF*R*h)/(J*var) %'
MF = inv(R*iXtX*R'); %'
mineig = min(eig(MF));
if(mineig <= 0)
  fprintf('ERROR: R*iXtX*Rt is not positive def (%g)\n',mineig);
  return;
end
CF = chol(MF);

% Compute Ct which makes the calucation of the t-stat easier:
% t = (Ct*R*h)/sqrt(var)
Mt = inv(diag(diag(R*iXtX*R'))); %'
mineig = min(eig(Mt));
if(mineig <= 0)
  fprintf('ERROR: R*diag(iXtX)*Rt is not positive def (%g)\n',mineig);
  return;
end
Ct = chol(Mt); % This is just a sqrt() because Mt is diagonal

% New DOF %
DOF = nruns * ntrs - size(X,2);

[nrows ncols nhstem fs ns endian bext] = fmri_bfiledim(hstem);
nvslice = nrows*ncols;
if(firstslice < 0) firstslice = fs; end
if(nslices < 0)    nslices = ns; end

lastslice = firstslice + nslices - 1;

%--------- Slice Loop ----------------------%
for nthslice = firstslice:lastslice
  fprintf('%2d ',nthslice);

  havg = [];
  hfile = sprintf('%s_%03d.%s',hstem,nthslice,bext);
  [havg eresvar sxadat] = fast_ldsxabfile(hfile);
  if(isempty(havg))
    fprintf('ERROR: could not load %s\n',hfile);
    return;
  end
  havg0 = havg;
  nhavg = size(havg,3);

  if(size(R,2) ~= nhavg)
    fprintf('ERROR: size(R,2) ~= nhavg (%d,%d)\n',size(R,2),nhavg);
    return;
  end

  havg = reshape(havg,[nvslice nhavg])'; %' nhavg X nvslice
  eresvar = reshape(eresvar,[nvslice 1])'; %' nhavg X nvslice

  % Adjust eresvar for new DOF %
  % eresvar = eresvar*hdr.DOF/DOF;

  %------------- FTest -------------------%
  if(~isempty(fstem) | ~isempty(fsigstem))
     Rhavg = R*havg;
     MeanRhavg = mean(Rhavg);
     SignMeanRhavg = sign(MeanRhavg);
     J = size(R,1);
     F = sum((CF*Rhavg).^2)./(J*eresvar);

     if(~isempty(fstem))
       tmp = reshape(F',[nrows ncols]); %'
       fname = sprintf('%s_%03d.bfloat',fstem,nthslice);
       fmri_svbfile(tmp,fname);
     end

     if(~isempty(fsigstem))
       fsig = FTest(J,DOF,F);
       ind = find(fsig < 10^(-80));
       fsig(ind) = 10^(-80);
       fsig = SignMeanRhavg .* -log10(-abs(fsig));
       tmp = reshape(fsig',[nrows ncols]); %'
       fname = sprintf('%s_%03d.bfloat',fsigstem,nthslice);
       fmri_svbfile(tmp,fname);
     end
  end % FTest %

  %------------- tTest -------------------%
  if(~isempty(tstem) | ~isempty(tsigstem) | ~isempty(tminsigstem))
     CtRhavg = (Ct*R)*havg;
     t = (CtRhavg)./repmat(sqrt(eresvar),[size(CtRhavg,1) 1]);

     if(~isempty(tstem))
       tmp = reshape(t',[nrows ncols size(t,1)]); %'
       fname = sprintf('%s_%03d.bfloat',tstem,nthslice);
       fmri_svbfile(tmp,fname);
     end

     if(~isempty(tsigstem) | ~isempty(tminsigstem))
       tsig = tTest(DOF,reshape1d(abs(t)),75);
       ind = find(tsig < 10^(-80));
       tsig(ind) = 10^(-80);
       tsig = sign(t) .* reshape(tsig,size(t));

       if(~isempty(tsigstem))
         tmp = reshape(tsig',[nrows ncols size(tsig,1)]); %'
         tmp = sign(tmp) .* -log10(abs(tmp));
         fname = sprintf('%s_%03d.bfloat',tsigstem,nthslice);
         fmri_svbfile(tmp,fname);
       end

       if(~isempty(tminsigstem)) 
         [tminsig iminsig] = min(abs(tsig));
         indminsig = sub2ind(size(tsig),iminsig,1:nvslice);
         tminsig = tsig(indminsig)*size(tsig,1); % Bonferroni
         tmp = reshape(tminsig',[nrows ncols]); %'
         tmp = sign(tmp) .* -log10(abs(tmp));
         fname = sprintf('%s_%03d.bfloat',tminsigstem,nthslice);
         fmri_svbfile(tmp,fname);
         % Could save imnisig here too
       end

     end
  end % tTest %

end % Slice Loop
fprintf('\n');

fprintf('stxgpredict finsished %g\n',toc);

