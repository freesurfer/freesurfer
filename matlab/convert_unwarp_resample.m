% $Id: convert_unwarp_resample.m,v 1.4 2003/08/12 18:14:25 ebeth Exp $
%
%% convert_unwarp_resample.m contains: 
% convert_unwarp_resample()
% load_dicom_plus()
% mdc()
% header2map(), type2map(), map2manuf() refer to TABLE = GRADWARPPATH/table.mat
%
%% convert_unwarp_resample() also calls:
% unwarp_init_globals() (GRADWARPPATH, TABLE, QuitOnError)
% unwarp_resample.m, load_mgh, save_mgh, etc.
%
%% also relevant:
% unwarp_scanners_table.m makes GRADWARPPATH/table.mat
%  - to change table.mat, add extra structure to unwarp_scanners_table.m and rerun
% $DEV/scripts/grad_unwarp
%   - invokes convert_unwarp_resample()
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_unwarp_resample(infile,series,outfile,corfovflag,unwarpflag,jacflag,interp_method,user_gradwarpfile,called_by_script)

% Converts dicom to mgh, possibly including unwarping and resampling
% to COR FOV.  Resampling to COR needs a rewrite.

% Can invoke with first 3 args = just save dicom as mgh
% For dewarping, need 5 or 6 or 7 or 8.

% Add stuff to path if necessary %
devdir = getenv('DEV');
% d = sprintf('%s/matlab',devdir);
% if(isempty(strfind(path,d))) path(path,d); end
d = sprintf('%s/fsfast/toolbox',devdir); % for qoe()
if(isempty(strfind(path,d))) path(path,d); end

global QuitOnError;
global GRADWARPPATH;
global TABLE;
if (~exist('called_by_script') | isempty(called_by_script))
  called_by_script=1; % default
end
unwarp_init_globals(called_by_script);

fprintf('GRADWARPPATH=%s\n',GRADWARPPATH);
fprintf('TABLE=%s\n',TABLE);
fprintf('QuitOnError=%d\n',QuitOnError);

tic;
% For if you're running matlab interactively instead of calling this
% from a script.

%%% check arguments %%%
if(nargin > 9 | nargin < 3)
  qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - 3-9 arguments required.');
  error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - 3-9 arguments required.');
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
  if strcmp(interp_method,'trilinear') interp_method ='linear'; end
  fprintf('Interp method is %s\n',interp_method);
end

if ((corfovflag | unwarpflag) & strcmp(interp_method,''))
  qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - If unwarping or resampling, must specify interpolation method');
  error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - If unwarping or resampling, must specify interpolation method');
end

%%% end check arguments %%%

%% load file, dicom or mgh %% 
if(isdicomfile(infile) | isdir(infile))
  fprintf('INFO: loading dicom volume at %s\n',infile);
  [vol, M0, mr_parms, Mdc, gradwarpfile] = load_dicom_plus(series,infile);
  if(isempty(vol))
    qoe('error loading dicom file.');
    error('error loading dicom file.');
  end
else % mgh
  fprintf('INFO: loading mgh-format infile %s\n',infile);
  % So now we know it's an mgh file - don't bother loading if args are wrong:
  if(unwarpflag & strcmp(user_gradwarpfile,''))
    qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - for unwarping an mgh file, user must supply gradwarpfilename or type');
    error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - for unwarping an mgh file, user must supply gradwarpfilename or type');
  end
  [vol, M0, mr_parms, Mdc] = load_mgh(infile);
  if(isempty(vol))
    qoe('error loading mgh file.');
    error('error loading mgh file.');
  end
end
%% end load file %% 

% That did or didn't give us a gradwarp file.  See if there's an override: % 
if(unwarpflag)
  if ~strcmp(user_gradwarpfile,'')
    if (isempty(strfind(user_gradwarpfile,'/')))
      % minor kludge: suggests user supplied a dewarp type (e.g. sonata)
      % instead of path (e.g. /space/dijon/foo/etc)
      gradwarpfile = type2map(user_gradwarpfile);
    else
      gradwarpfile = user_gradwarpfile;
    end
  elseif (~exist('gradwarpfile') | isempty(gradwarpfile) | strcmp(gradwarpfile,''))
    qoe('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - no gradwarpfile found - user must supply gradwarpfile or type');
    error('convert_unwarp_resample(infile,serno,outfile,corfovflag,unwarpflag,jacflag,interp_method,gradwarpfile) - no gradwarpfile found - user must supply gradwarpfile or type');
  end
  [placeholder,inflag,thruflag]=map2manuf(gradwarpfile);
end

% If there wasn't much to do:
if(~corfovflag & ~unwarpflag)
  fprintf('Writing MGH output file %s\n',outfile);
  save_mgh(vol,outfile,M0,mr_parms);
  return;
end

% Convert M to 1-based because that is what anders code expects
M = vox2ras_0to1(M0);

% Setup the output geometry %
if(corfovflag)
  % This might need to be modified so as to center the volume
  Mout = [0    -1     0   129;...
          0     0     1  -129;...
         -1     0     0   129;...
          0     0     0     1];
  size_out = [256 256 256];
else 
   Mout = M;
  size_out = size(vol);
end

% Do unwarping and/or resampling %
fprintf('Beginning unwarping and/or resampling \n');

[volout, Mout] = unwarp_resample(vol,M,size_out,Mout,Mdc,unwarpflag,jacflag,0,interp_method,inflag,thruflag,gradwarpfile);

% Convert Mout back to 0-based %
M0out = vox2ras_1to0(Mout);

fprintf('Writing MGH output file %s\n',outfile);
save_mgh(volout,outfile,M0out,mr_parms);

fprintf('convert_unwarp_resample done (%g)\n',toc);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol, M, mr_parms, Mdc, gradfilename] = load_dicom_plus(seriesno,dcm)

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

% Add stuff to path if necessary

global QuitOnError;

devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir);
if(isempty(findstr(d,path))) path(path,d); end

if(nargin < 1 | nargin > 3)
  qoe('[vol, M, mr_parms, Mdc, gradfilename] = load_dicom_plus(<seriesno or ''>, <dcmdir or dcmfile>)');
  error('[vol, M, mr_parms, Mdc, gradfilename] = load_dicom_plus(<seriesno or ''>, <dcmdir or dcmfile>)');
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

if isfield(dcminfo(1),'ScannerSerialNo')
  SSN = dcminfo(1).ScannerSerialNo;
else SSN = ''; end

if isfield(dcminfo(1),'InstitutionName')
  IN = dcminfo(1).InstitutionName;
else IN = ''; end

if isfield(dcminfo(1),'StationName')
  SN = dcminfo(1).StationName;
else SN = ''; end

gradfilename = header2map(dcminfo(1).Manufacturer,dcminfo(1).ManufacturersModelName,SSN,IN,SN);

Mdc = mdc(dcminfo);

return;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mdc = mdc(dcminfo)
% Matrix of direction cosines

Mdc = zeros(3,3);
Mdc(:,1) = dcminfo(1).ImageOrientationPatient(1:3);
Mdc(:,2) = dcminfo(1).ImageOrientationPatient(4:6);
sdc = dcminfo(2).ImagePositionPatient-dcminfo(1).ImagePositionPatient;
Mdc(:,3) = sdc / sqrt(sum(sdc.^2));
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >> header2map('siemens','sonata','21913','','')

function gradfilename = header2map(M,MMN,SSN,IN,SN)

global TABLE
load(TABLE);

INSN = sprintf('%s %s',IN,SN);

for ii=[1:length(m)]
  if (strcmpi(m(ii).Manufacturer,M) & any(strcmpi(m(ii).ManufacturersModelName,MMN)))
    if(strcmpi(m(ii).Manufacturer,'siemens'))
      % For Siemens, M and MMN fully determine the gradfilename
      gradfilename = m(ii).filename;
      return;
    end
    if (any(strcmpi(m(ii).ScannerSerialNo,SSN)) | any(strcmpi(m(ii).InstitutionStation,INSN)))
      % For GE, need either of these to narrow it down.
      gradfilename = m(ii).filename;
      return;
    end
  end
end

gradfilename = '';
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradfilename = type2map(type)

global TABLE
load(TABLE);

% This didnt work:
% map = m(find(strcmp(m.nickname,type)));
% gradfilename= map(1);

for ii=[1:length(m)]
  if (strcmp(m(ii).nickname,type))
    gradfilename = m(ii).filename;
    return
  end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [manuf, inflag, thruflag] = map2manuf(gradfilename)

global TABLE
load(TABLE);

manuf = ''; inflag = []; thruflag = [];
gfn = gradfilename(max(strfind(gradfilename,'/'))+1:end);  % i.e. after last "/"

for ii=[1:length(m)]
  testgfn = m(ii).filename(max(strfind(gradfilename,'/'))+1:end);
  if (strcmpi(testgfn,gfn))
    manuf = lower(m(ii).Manufacturer);
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

