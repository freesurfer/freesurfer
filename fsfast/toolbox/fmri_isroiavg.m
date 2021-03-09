% fmri_isroiavg.m
%   Takes the output of fmri_roiavg for each subject and computes
%   inter-subject statistics using a fixed effects model
%
%
%
% DelaySign, Contrast
%


%
% fmri_isroiavg.m
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

fprintf('\n$Id: fmri_isroiavg.m,v 1.3 2011/03/02 00:04:06 nicks Exp $\n');
nROIFiles = length(roifile);

FilesUsed = [];

hAllSubj = [];
for n = 1:nROIFiles
  infile = deblank(roifile(n).name);
  %fprintf('Loading %s\n',infile);
  load(infile);
  %fprintf('%s %d\n',infile,dof(1));
  if(dof(1) ~= 0)
    hSubjAvg = reshape1d(hAvg');
    hAllSubj = [hAllSubj hSubjAvg];
    FilesUsed = strvcat(FilesUsed,infile);
  else
    fprintf('File %s has no samples\n',infile);
  end
end

Nh = length(DelaySign);
Nc = length(Contrast);

[nD nSubj] = size(hAllSubj);
nD2 = Nh*Nc;
if(nD ~= nD2)
  msg = sprintf(...
   'ERROR: DelaySign or Contrast does not have correct number of components (%d/%d)',...
                nD,nD2);
  qoe(msg);error(msg);
end

v = reshape1d( DelaySign' * Contrast );
if(DSAvg)
  mm = 0;
  ind = find(v == 1);
  if(~isempty(ind)) 
     v(ind) = v(ind)/length(ind); 
     mm = mm + 1;
  end

  ind = find(v == -1);
  if(~isempty(ind)) 
     v(ind) = v(ind)/length(ind); 
     mm = mm + 1;
  end
  v = v/mm;
else
  ind = find(v ~= 0);
  v = v/length(ind);
end

v2 = repmat(v, [1 nSubj]);
s = sum(v2 .* hAllSubj);   % this is actually the mean %

dof = nSubj - 1;
isavg = mean(s);
isstd = std(s);
iserr = isstd/sqrt(dof);
T = isavg/iserr;
p = sign(T)*tTest(dof,T);

fprintf('Contrast: ');
fprintf('%d ',Contrast);
fprintf('\n');
fprintf('DelaySign:');
fprintf('%d ',DelaySign);
fprintf('\n');
for n = 1:nSubj
  fprintf('%2d %s  %10.4f\n',n,FilesUsed(n,:),s(n));
end

fprintf('\n area dof = %d, avg = %g, stdev = %g, stderr = %g, t = %g, p = %g\n\n',...
        dof, isavg,isstd,iserr,T,p);

if(report)
  fid = fopen(reportfile,'w');
  fprintf(fid,'------------------------------------------\n');
  fprintf(fid,'Threshold = %g\n',threshold);
  fprintf(fid,'ROI (slc,row,col) = (%d,%d), (%d,%d), (%d,%d)\n',...
               FirstSlice,FirstSlice+nSlices-1,roi_rowcol(1),roi_rowcol(3),...
               roi_rowcol(2),roi_rowcol(4));
  fprintf(fid,'Contrast: ');
  fprintf(fid,'%d ',Contrast);
  fprintf(fid,'\n');
  fprintf(fid,'DelaySign:');
  fprintf(fid,'%d ',DelaySign);
  fprintf(fid,'\n');
  for n = 1:nSubj
    fprintf(fid,'%2d %s  %10.4f\n',n,FilesUsed(n,:),s(n));
  end

  fprintf(fid,'\ndof = %d\navg = %g\nstdev = %g\nstderr = %g\nt = %g\np = %g\n\n',...
        dof, isavg,isstd,iserr,T,p);
  fclose(fid);
end

if(ShowResults)
  hAllAvg = mean(hAllSubj');
  hAllStd = std(hAllSubj');
  hAllAvg2 = reshape(hAllAvg,[Nh Nc]);
  hAllStd2 = reshape(hAllStd,[Nh Nc]);
  dof1  = dof*ones(hdrdat.Nc,1);
  hHDR = figure(1);
  hdrviewlite('init',hHDR,t,dof1);
  hdrviewlite('plot',hHDR,hAllAvg2',hAllStd2', [0 0]);
  if(ShowResults ~= 2)  uiwait(hHDR); end
end

return

