function dof = fmri_lddof(varargin)
%
% dof = fmri_lddof(BFileName)
%
% Loads info from the .dof file associated with BFileName.
% This is only relevant for back-compatibility with selavg.
%
%


%
% fmri_lddof.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

dof = [];

if(nargin == 0) 
  fprintf(2,'USAGE: fmri_lddof(BFileName)');
  qoe;
  return;
end

if( length(varargin) == 1)
  BFileList = varargin{1};
  nRuns = size(BFileList,1);
else
  nRuns = length(varargin);
  BFileList = '';
  for r = 1:nRuns,
    BFileList = strvcat(BFileList,varargin{r});
  end
end


for r = 1:nRuns,

  BFileName = deblank(BFileList(r,:));
  ks = findstr(BFileName,'.bshort');
  kf = findstr(BFileName,'.bfloat');

  if(isempty(ks) & isempty(kf))
    msg = 'BFileName must be either bshort or bfloat';
    qoe(msg);
    error(msg);
  end

  if( ~isempty(ks) ) 
    Base = BFileName(1:ks-1);
  else               
    Base = BFileName(1:kf-1);
  end

  if( isempty(Base) )
    s = 'BFileName must have a non-null base';
    qoe(msg);
    error(msg);
  end

  %%% Open the header file %%%%
  DOFFile = strcat(Base,'.dof');
  fid=fopen(DOFFile,'r');
  if fid == -1 
    msg = sprintf('Could not open %s file',DOFFile); 
    qoe(msg);
    error(msg);
  end

  %%%% Read the Dimension from the header %%%%
  d=fscanf(fid,'%g');
  fclose(fid);

  n = length(d);
  nCond = n/3;
  d = reshape(d,[3 nCond])';
  dof = d(:,[2 3]);

end

return;

%%% y now has size(y) = [nR nC nD nRuns] %%%

