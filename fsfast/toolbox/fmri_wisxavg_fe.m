% Name: fmri_isxavg_fe
% Purpose: implements fixed-effects intersession averaging
%          for output of selxavg
% Author: Douglas Greve
% Questions or Comments: analysis-bugs@nmr.mgh.harvard.edu

%%%% These variables must be specified %%%%%%%%
% InStemList
% FirstSlice
% NSlices
% OutStem


%
% fmri_wisxavg_fe.m
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

if( exist(deblank('InStemList')) ~= 1)
  msg = sprintf('Variable InStemList does not exist');
  qoe(msg);error(msg);
end

nsessions = size(InStemList,1);
LastSlice = FirstSlice + NSlices - 1;

for slice = FirstSlice:LastSlice

  % fprintf('Processing Slice %d\n',slice);
  fprintf('%2d ',slice);
  eVarSum = 0;
  DOFSum  = 0;
  hAvgSum = 0;
  hCovSum = 0;

  for session = 1:nsessions,
    % fprintf('Session %d\n',session);
    InStem  = deblank(InStemList(session,:));
    InSA    = sprintf('%s_%03d.bfloat',InStem,slice);
    DatFile = sprintf('%s.dat',InStem);
    InHOffset  = sprintf('%s-offset_%03d.bfloat',InStem,slice);
    
    hd = fmri_lddat3(DatFile);
    ysa  = fmri_ldbfile(InSA);

    if(pctsigch)
      hoffset = fmri_ldbfile(InHOffset);
      hofftmp = repmat(hoffset,[1 1 size(ysa,3)]);
      ysa = ysa./hofftmp;
    end

    [hAvg eVar] = fmri_sa2sxa(ysa,hd.Nh);
    hCov = hd.hCovMtx;
    eVarSum = eVarSum + eVar * hd.DOF;
    hAvgSum = hAvgSum + hAvg * hd.DOF;
    hCovSum = hCovSum + inv(hCov);
    DOFSum  = DOFSum + hd.DOF;
  end % loop over sessions %

  hAvgGrp = hAvgSum/DOFSum;
  eVarGrp = eVarSum/DOFSum;
  hCovGrp = inv(hCovSum);

  hd.hCovMtx = hCovGrp;
  hd.DOF = DOFSum;
  [ySA dof] = fmri_sxa2sa(eVarGrp,hCovGrp,hAvgGrp,hd);

  OutSA = sprintf('%s_%03d.bfloat',OutStem,slice); 
  fmri_svbfile(ySA, OutSA); 

  OutDat = sprintf('%s.dat',OutStem);
  fmri_svdat2(OutDat,hd);

%  OutCov = sprintf('%s_hcov_%03d.bfloat',OutStem,slice); 
%  fmri_svbfile(hCovGrp,OutCov);

%%%% Dont save dof file anymore %%%%%%%%%%%%%%%%%%
%  dofFile = sprintf('%s_%03d.dof',OutStem,slice); 
%  fid=fopen(deblank(dofFile),'w');
%  if( fid == -1 )
%    msg = sprintf('Could not open dof file %s\n',dofFile);
%    qoe(msg);  error(msg);
%  end
%  for c = 0:hd.Nnnc,
%    fprintf(fid,'%5d %5d %5d\n',c, dof(c+1), dof(c+1)-1);    
%  end
%  fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % Loop over slices %
fprintf('\n');

fprintf(1,'fmri_isavg_fe completed SUCCESSUFLLY\n\n');
