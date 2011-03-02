function convert_unwarp_resample(infile,series,outfile,corfovflag,unwarpflag,jacflag,interp_method,user_gradwarpfile,called_by_script)
%
%% convert_unwarp_resample.m contains: 
% convert_unwarp_resample()
% load_dicom_and_stuff()
% mdc()
% header2map(), type2map(), map2manuf() refer to TABLE = GRADWARPPATH/table.mat
%
%% convert_unwarp_resample() also calls:
% unwarp_init_globals.m (initializes GRADWARPPATH, TABLE, QuitOnError)
% unwarp_resample.m, load_mgh2.m, save_mgh.m, load_dicom_series.m, etc.
%
%% also relevant:
% unwarp_scanners_table.m creates GRADWARPPATH/table.mat
%  - to change table.mat, edit unwarp_scanners_table.m structure and rerun it
% $DEV/scripts/grad_unwarp
%   - invokes convert_unwarp_resample()
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% convert_unwarp_resample.m
%
% Original Author: Elizabeth Haley
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.13 $
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



% Converts dicom to mgh, possibly including unwarping and resampling
% to COR FOV.  Resampling to COR needs a rewrite.

% Can invoke with first 3 args = just save dicom as mgh
% For dewarping, need 5 or 6 or 7 or 8.

% Add stuff to path if necessary %
devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir); % for qoe()
if(isempty(strfind(path,d))) path(path,d); end

global GRADWARPPATH;
global TABLE;
global QuitOnError;
% QuitOnError is for: if you're running matlab interactively instead
% of calling this from a script and there's an error, it fails
% gracefully (matlab quits).  unwarp_init_globals sets global
% QuitOnError based on value of called_by_script.
% convert_unwarp_resample.m is usually called by script grad_unwarp.
if (~exist('called_by_script') | isempty(called_by_script))
  called_by_script=1; % default
end
unwarp_init_globals(called_by_script);

fprintf('GRADWARPPATH=%s\n',GRADWARPPATH);
fprintf('TABLE=%s\n',TABLE);
fprintf('QuitOnError=%d\n',QuitOnError);

tic;

%%% check arguments %%%
if(nargin > 9 | nargin < 3)
  qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - 3-9 arguments required.');
  error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - 3-9 arguments required.');
end

if (~exist('corfovflag') | isempty(corfovflag)) corfovflag=0; end
if (~exist('unwarpflag') | isempty(unwarpflag)) unwarpflag=0; end
if (~exist('jacflag')    | isempty(jacflag))    jacflag=1;    end
if (~exist('user_gradwarpfile') | isempty(user_gradwarpfile))
  user_gradwarpfile='';
end
if (~exist('interp_method') | isempty(interp_method))
  interp_method='';
else
  % Other names for 'linear' sampling:
  if (strcmp(interp_method,'trilinear')|strcmp(interp_method,'interpolate')) interp_method ='linear'; end
  fprintf('Interp method is %s\n',interp_method);
end

if ((corfovflag | unwarpflag) & strcmp(interp_method,''))
  qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - If unwarping or resampling, must specify interpolation method');
  error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - If unwarping or resampling, must specify interpolation method');
end
%%% end check arguments %%%

%% load file, dicom or mgh %% 
if((exist(infile,'file')) & ...
   ((strcmp(infile(end-3:end),'.mgh')) | (strcmp(infile(end-3:end),'.mgz'))))
  % Mildly cheesy test for mgh file, good enough.
  fprintf('INFO: loading mgh-format infile %s\n',infile);
  % So now we know it's an mgh file - dewarping an mgh file requires
  % user-supplied gradwarpfile or gradwarptype - don't bother loading
  % if args are wrong:
  if(unwarpflag & strcmp(user_gradwarpfile,''))
    qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - for unwarping an mgh file, user must supply gradwarpfilename or type');
    error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - for unwarping an mgh file, user must supply gradwarpfilename or type');
  end
  fprintf('INFO: loading mgh volume %s\n',infile);
  [vol, M0, mr_parms] = load_mgh(infile);
  if(isempty(vol))
    qoe('error loading mgh file.');
    error('error loading mgh file.');
  end
elseif (isdicomfile(infile) | isdir(infile))
  % dicom file or dicom directory
  fprintf('INFO: loading dicom volume at %s\n',infile);
  [vol, M0, mr_parms, gradwarpfile] = load_dicom_and_stuff(series,infile);
  if(isempty(vol))
    qoe('error loading dicom file.');
    error('error loading dicom file.');
  end
else
  qoe('convert_unwarp_resample: infile was neither dicom nor mgh?');
  error('convert_unwarp_resample: infile was neither dicom nor mgh?');
end
%% end load file %% 

% That did or didn't give us a gradwarp file.  See if there's an override: % 
if(unwarpflag)
  if ~strcmp(user_gradwarpfile,'')
    % Minor kludge: suggests user supplied a dewarp type (e.g. sonata)
    % instead of path (e.g. /space/dijon/foo/etc).  If it's not a
    % dewarp type, type2map will notice
    if (isempty(strfind(user_gradwarpfile,'/')))
      gradwarpfile = type2map(user_gradwarpfile);
    else
      gradwarpfile = user_gradwarpfile;
    end
  elseif (~exist('gradwarpfile') | isempty(gradwarpfile) | strcmp(gradwarpfile,''))
    qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - no gradwarpfile found - user must supply gradwarpfile or type');
    error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile,called_by_script) - no gradwarpfile found - user must supply gradwarpfile or type');
  end

  % If unwarping, have to figure out whether inplane or thruplane or
  % not (no matter where the gradwarpfile came from (user-supplied
  % filename, user-supplied type, dicom headers)) - so look the
  % filename up in table.mat.  Match gradwarpfile to the right data
  % struct, find the manufacturer, and, from the manufacturer,
  % conclude whether volume has already been inplane-dewarped (like
  % GE).
  [manuf_dontcare,inflag,thruflag]=map2manuf(gradwarpfile);
else
  inflag=0; thruflag=0; gradwarpfile='';
end

% If there wasn't much to do (e.g. just converting dicom to mgh, no
% unwarp or resample):
if(~corfovflag & ~unwarpflag)
  fprintf('Writing MGH output file %s\n',outfile);
  save_mgh(vol,outfile,M0,mr_parms);
  return;
end

% We'll need ras_xform xyzc_ras, Mdc for GE throughplane projection
% and c_ras for centering cor volume, if corfovflag:
[Mdc, c_ras] = ras_xform_from_vox2ras_xform(M0,size(vol));

% Convert M to 1-based because that is what anders code expects
M = vox2ras_0to1(M0);

% Set up the output geometry %
if(corfovflag)
  %%%%%
  % Q: This might need to be modified so as to center the volume
  % A: Yeah, this centers it on c_ras = 0,0,0 - scanner isocenter!
  %%%%%
  % Q: Maybe before I messed up Anders's code, this matrix ensured
  % cor orientation - but at this point it surely does not.  However,
  % if you want COR format, you're probably mri_converting anyway, and
  % you can start with this matrix and mri_convert will take care of
  % it with no interpolation, no problem.
  %   Mout = [0    -1     0   129;...
  %           0     0     1  -129;...
  %          -1     0     0   129;...
  %           0     0     0     1];
  %%%%%
  
  Mout0 = [-1     0     0   128+c_ras(1) ;...
            0     0     1  -128+c_ras(2) ;...
            0    -1     0   128+c_ras(3) ;...
            0     0     0     1];
  %% WORKS!
  % Convert Mout to 1-based because that is what anders code expects
  Mout = vox2ras_0to1(Mout0);
  size_out = [256 256 256];
else 
  Mout = M; % Already converted from M0 to 1-based
  size_out = size(vol);
end

% Do unwarping and/or resampling %
fprintf('Beginning unwarping and/or resampling \n');

% but unwarp_resample does not alter Mout
[volout, Mout] = unwarp_resample(vol,M,size_out,Mout,Mdc,unwarpflag,jacflag,0,interp_method,inflag,thruflag,gradwarpfile);

% Convert Mout back to 0-based %
M0out = vox2ras_1to0(Mout);

fprintf('Writing MGH output file %s\n',outfile);
save_mgh(volout,outfile,M0out,mr_parms);

fprintf('convert_unwarp_resample done (%g)\n',toc);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [vol, M, mr_parms, Mdc, c_ras, gradfilename] = load_dicom_and_stuff(seriesno,dcm)
function [vol, M, mr_parms, gradfilename] = load_dicom_and_stuff(seriesno,dcm)

% Reads in a dicom series given:
%  1. The series number and directory, or
%  2. nothing and a dicom file from the desired series
%
% If the series number is given but no dcmdir is given, then the
% current directory is assumed. All files in the dcmdir are examined
% and the dicom files for the given series are then loaded.
%
% If a dicom file is given, then seriesno and dcmdir are determined
% from the file and file name.
%
% mr_parms = [tr flipangle te ti]
%
% Bugs: will not load multiple frames or mosaics properly.

global QuitOnError;

% Add stuff to path if necessary
devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir); % for qoe
if(isempty(findstr(d,path))) path(path,d); end

if(nargin < 1 | nargin > 3)
  qoe('[vol, M, mr_parms, Mdc, gradfilename] = load_dicom_and_stuff(<seriesno or ''>, <dcmdir or dcmfile>)');
  error('[vol, M, mr_parms, Mdc, gradfilename] = load_dicom_and_stuff(<seriesno or ''>, <dcmdir or dcmfile>)');
end

if (isdicomfile(dcm))
  [vol M dcminfo mr_parms] = load_dicom_series('','',dcm); % dcm was a file
elseif isdir(dcm)
  [vol M dcminfo mr_parms] = load_dicom_series(seriesno,dcm); % dcm was a dir
else
  qoe('dcm was neither a dicom file nor a directory.');
  error('dcm was neither a dicom file nor a directory.');
end;

if(isempty(vol)) return; end

% Figure out what the appropriate gradwarp file is for this machine,
% based on ScannerSerialNo or DeviceSerialNumber, InstitutionName, and
% StationName.

if (isfield(dcminfo(1),'ScannerSerialNo') & ~isempty(dcminfo(1).ScannerSerialNo))
  SSN = dcminfo(1).ScannerSerialNo;
elseif (isfield(dcminfo(1),'DeviceSerialNumber') & ~isempty(dcminfo(1).DeviceSerialNumber))
  SSN = dcminfo(1).DeviceSerialNumber;
else SSN = ''; end

if isfield(dcminfo(1),'InstitutionName')
  IN = dcminfo(1).InstitutionName;
else IN = ''; end

if isfield(dcminfo(1),'StationName')
  SN = dcminfo(1).StationName;
else SN = ''; end

gradfilename = header2map(dcminfo(1).Manufacturer,dcminfo(1).ManufacturerModelName,SSN,IN,SN);
fprintf('load_dicom_and_stuff: INFO: gradwarpfile is %s\n',gradfilename); %EDEBUG%

return;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >> header2map('siemens','sonata','21913','','')

function gradfilename = header2map(M,MMN,SSN,IN,SN)

global GRADWARPPATH;
global TABLE;
load(TABLE);

INSN = sprintf('%s %s',IN,SN);

gradfilename='';
for ii=[1:length(m)]
  if (strcmpi(m(ii).Manufacturer,M) & any(strcmpi(m(ii).ManufacturersModelName,MMN)))
    if(strcmpi(m(ii).Manufacturer,'siemens'))
      % For Siemens, M and MMN fully determine the gradfilename.  So far.  We think.
      % gradfilename = m(ii).filename;
      gradfilename = sprintf('%s/%s',GRADWARPPATH,m(ii).filename);
      return;
    end
    if (any(strcmpi(m(ii).ScannerSerialNo,SSN)) | any(strcmpi(m(ii).InstitutionStation,INSN)))
      % For GE, need either of these to narrow it down.
      % gradfilename = m(ii).filename;
      gradfilename = sprintf('%s/%s',GRADWARPPATH,m(ii).filename);
      return;
    end
  end
end

gradfilename = '';
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradfilename = type2map(type)

global GRADWARPPATH;
global TABLE;
load(TABLE);

% This didnt work:
% map = m(find(strcmp(m.nickname,type)));
% gradfilename= map(1);

for ii=[1:length(m)]
  if (strcmp(m(ii).nickname,type))
    % gradfilename = m(ii).filename;
    gradfilename = sprintf('%s/%s',GRADWARPPATH,m(ii).filename);
    return
  end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [manuf, inflag, thruflag] = map2manuf(gradfilename)

global GRADWARPPATH;
global TABLE;
load(TABLE);

manuf = ''; inflag = []; thruflag = [];
gfn = gradfilename(max(strfind(gradfilename,'/'))+1:end);  % i.e. after last "/"

for ii=[1:length(m)]
  % Changed this: m.filenames now have no path.
  % testgfn = m(ii).filename(max(strfind(gradfilename,'/'))+1:end);
  testgfn = m(ii).filename;
  % Find the entry with the same gradfilename %
  if (strcmpi(testgfn,gfn))
    manuf = lower(m(ii).Manufacturer);
    % Determine dewarp projection from manufacturer - assume GE is
    % always already in-plane dewarped and that Siemens is not.
    if strcmpi(manuf,'siemens')
      inflag=0; thruflag=0; % volume to be full-3D dewarped.
    elseif strcmpi(manuf,'ge medical systems');
      inflag=0; thruflag=1; % volume to be throughplane-only dewarped.
    else
      fprintf('INFO: dewarping dimensions (throughplane/full3D) not set.\n');
    end
    return;
  end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mdc, c_ras] = ras_xform_from_vox2ras_xform(M,size_vol)

MdcD = M(1:3,1:3);
delta = sqrt(sum(MdcD.^2));
Mdc = MdcD./repmat(delta,[3 1]);
Pcrs_c = [size_vol(1)/2; size_vol(2)/2; size_vol(3)/2; 1];
Pxyz_c = M*Pcrs_c;
c_ras = Pxyz_c(1:3);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

