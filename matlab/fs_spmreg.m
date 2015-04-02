function err = fs_spmreg(targvol,srcvol,outvol,DOFStr,costfun,errfile,UseSPMGetSpaceStr,fmt,monlyStr)
% err = fs_spmreg(targvol,srcvol,outvol,DOFStr,errfile,UseSPMGetSpaceStr,fmt,monlyStr)
% matlab script to run spm registration. It is set up to be %
% "deployable" so all the args must be strings.
% 1 targvol - name of target/reference volume
% 2 srcvol - name of moveable/source volume
% 3 outvol - name of output volume (with UseSPMGetSpace=0)
% 4 DOFStr - registration DOF: 6, 9, 12
% 5 costfun - cost function (nmi,mi,ecc,ncc)
% 6 errfile - file created when there is an error; the creation of
%   this file is really the only way to know if there was an error
% 7 UseSPMGetSpaceStr - 0=dont, 1 = use MRIwrite()
% 8 fmt - format (nii, nii.gz, or img) only used for error checking
% 9 monlyStr - only controls whether errfile is created on error (0=created) 

err = 1;
fprintf('Staring fs_spmreg\n');

DOF = sscanf(DOFStr,'%d',1);
UseSPMGetSpace = sscanf(UseSPMGetSpaceStr,'%d',1);
monly = sscanf(monlyStr,'%d',1);
fprintf('using %s as the target image\n',targvol);
fprintf('using %s as the source image\n',srcvol);
fprintf('output is %s \n',outvol);
fprintf('DOF %d\n',DOF);
fprintf('costfun %s\n',costfun);
fprintf('UseSPMGetSpace %d\n',UseSPMGetSpace);
fprintf('fmt %s\n',fmt);
fprintf('errfile %s\n',errfile);

if(exist('spm_coreg') ~= 2)
  fprintf('ERROR: cannot find spm_coreg.m. Make sure that the SPM\n');
  fprintf('   package is in your matlab path (check ~/matlab/startup)\n');
  fp = fopen(errfile,'w');
  fprintf(fp,'ERROR: cannot find spm_coreg.m. Make sure that the SPM\n');
  fprintf(fp,'   package is in your matlab path (check ~/matlab/startup)\n');
  fclose(fp);
  return; 
end

which spm_coreg

% This is a hack to try to determine which version of spm is being run. If
% spm8 is being run with analyze as input, a left-right reveral will occur.
% This searches for "spm8" in the spm path, not perfect, but hopefully good 
% enough. It does not appear that spm offers a "version" command.
IsSPM8 = length(findstr(dirname(which('spm_coreg')),'spm8'));
if(IsSPM8 & strcmp(fmt,'img'))
  fprintf('\n\n');
  fprintf('ERROR: you appear to be using spm8. If so, re-run this with --nii\n');
  fprintf('\n\n');
  fp = fopen(errfile,'a');
  fprintf(fp,'ERROR: you appear to be using spm8. If so, re-run this with --nii\n');
  fclose(fp);
  return; 
end

global defaults
spm_defaults % load spm2's default fields

fprintf('INFO: Assuming RAS coordinate system\n');
defaults.analyze.flip = 0;

%===========================================================================
% coregistration  defaults
%===========================================================================
params = [0 0 0 0 0 0];
if(DOF ==  9) params = [params 1 1 1]; end
if(DOF == 12) params = [params 1 1 1 0 0 0]; end

defaults.coreg.estimate.cost_fun = costfun;
defaults.coreg.estimate.sep      = [4 2];
defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
defaults.coreg.estimate.fwhm     = [7 7];
defaults.coreg.estimate.params	 = [params];
defaults.coreg.estimate.graphics = 0; % dont create ps file

%===========================================================================
% coregister files
%===========================================================================


% This can fail if there is not a path for each volume (ie, include
% ./ if current folder)
[pth1,nam1,ext1] = fileparts(deblank(targvol));
[pth2,nam2,ext2] = fileparts(deblank(srcvol));
tmp = sprintf('%s/%s%s',pth1,nam1,ext1);
%tmp = spm_get('Files',pth1,[nam1 ext1]);
VG = spm_vol(tmp);	% target image
if(isempty(VG))
  fprintf('ERROR: loading target %s\n',targvol);
  fp = fopen(errfile,'a');
  fprintf(fp,'ERROR: loading target %s\n',targvol);
  fclose(fp);
  return;
end
% tmp = spm_get('Files',pth2,[nam2 ext2])
tmp = sprintf('%s/%s%s',pth2,nam2,ext2);
VF = spm_vol(tmp);	% source image
if(isempty(VF))
  fprintf('ERROR: loading source %s\n',srcvol);
  fp = fopen(errfile,'a');
  fprintf(fp,'ERROR: loading source %s\n',srcvol);
  fclose(fp);
  return; quit;
end

fprintf('determining coregistration parameters...\n');
fprintf('\n\nINFO: ignore warnings about scalefactor\n\n');

tic;
try
  x  = spm_coreg(VG.fname,VF.fname,defaults.coreg.estimate);
catch ME
  fprintf('ERROR: spm_coreg\n');
  ME
  fp = fopen(errfile,'a');
  fprintf(fp,'ERROR: spm_coreg (catch)\n');
  fprintf(fp,'%s\n',which('spm_coreg'));
  fprintf(fp,'%s\n',ME.message);
  fclose(fp);
  return;
end
if(isempty(x))
  fprintf('ERROR: spm_coreg\n');
  if(~monly) 
    fp = fopen(errfile,'a');
    fprintf(fp,'ERROR: spm_coreg, x is empty\n');
    fclose(fp);
    return; quit; 
  end
end
fprintf('Finished after %g sec\n',toc);

fprintf('Parameters ');
fprintf('%g ',x);
fprintf('\n');

%---------write out the .mat file---see spm_coreg_ui.m line 283

M  = inv(spm_matrix(x));

if(UseSPMGetSpace)
  fprintf('Using spm_get_space\n');
  MM = zeros(4,4,1);
  MM = spm_get_space(deblank(VF.fname));
  spm_get_space(deblank(VF.fname), M*MM);
  fprintf('\n\nINFO: ignore warnings about scalefactor\n\n');
else
  fprintf('Using MRIread/write instead of spm_get_space\n');
  ff = MRIread(deblank(VF.fname));
  ff.vox2ras1 = M*ff.vox2ras1;
  ff.vox2ras0 = vox2ras_1to0(ff.vox2ras1);
  ff.vox2ras = ff.vox2ras0;
  MRIwrite(ff,outvol);
end

err = 0;
fprintf('fs_spmreg done\n');

return;
