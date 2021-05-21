function cmtx = fast_hypospec2cmtx(hsfile)
%
% cmtx = fast_hypospec2cmtx(hsfile)
%
% Converts a hypothesis specification file into
% a contrast matrix.


%
% fast_hypospec2cmtx.m
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

cmtx = [];

if(nargin ~= 1)
  msg = 'USAGE: cmtx = fast_hypospec2cmtx(hsfile)';
  qoe(msg); error(msg);
end

fid = fopen(hsfile);
if(fid == -1)
  msg = sprintf('Error opening %s',hsfile);
  qoe(msg); error(msg);
end

nvariates = readkeyvalue(fid,'nvariates',1,'%d');
if(isempty(nvariates))
  fprintf('ERROR reading %s\n',hsfile);
  fprintf('Could not read nvariates\n');
  return;
end

nconditions = readkeyvalue(fid,'nconditions',1,'%d');
if(isempty(nconditions))
  fprintf('ERROR reading %s\n',hsfile);
  fprintf('Could not read nconditions\n');
  return;
end

nestspercond = readkeyvalue(fid,'nestspercond',nconditions,'%d');
if(isempty(nestspercond))
  fprintf('ERROR reading %s\n',hsfile);
  fprintf('Could not read nestspercond\n');
  return;
end

condweights = readkeyvalue(fid,'conditionweights',nconditions,'%f');
if(isempty(condweights))
  fprintf('ERROR reading %s\n',hsfile);
  fprintf('Could not read condweights\n');
  return;
end

nthcondspec = 1;
while(~feof(fid))
  
  tmp = fgetnextline(fid);
  if(isempty(tmp)) break; end
  tmp = sscanf(tmp,'%s',1);
  if(strcmp('tmp','condspec'))
    fprintf('ERROR reading condspec %s from %s\n',nthcondspec,hsfile);
    return;
  end

  [CSM CSL] = ParseCondSpec(fid,nvariates,nconditions,nthcondspec);
  if(isempty(CSM)) 
    fprintf('ERROR reading condspec %s from %s\n',nthcondspec,hsfile);
    return; 
  end

  CondSpec(nthcondspec).Mtx = CSM;
  CondSpec(nthcondspec).List = CSL;
  nthcondspec = nthcondspec + 1;

end % while(~feof(fid))

ncondspecs = nthcondspec - 1;
if(ncondspecs == 0)
  fprintf('ERROR: no condition specs in %s\n',hsfile);
  return;
end

cm = [];
for cond = 1:nconditions
  if(condweights(cond) ~= 0)
    i = [];
    for cs = 1:ncondspecs,
      CSL = CondSpec(cs).List;
      ii = find(CSL == cond);
      if(~isempty(ii)) i = [i cs]; end
      %fprintf('%d  %d\n',cond,ii);
    end
    if(isempty(i))
      fprintf('ERROR: condition %d does not have a condspec\n',cond);
      return;
    end
    if(length(i) > 1)
      fprintf('ERROR: condition %d has multiple condspecs\n',cond);
      return;
    end
    m = CondSpec(i).Mtx;
  else
    m = zeros(nvariates,nestspercond(cond));
  end

  cm = [cm condweights(cond)*m];
end

cmtx = cm;

return
%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
%----------------------------------------------------------------------%


%----------------------------------------------------------------------%
function [CondSpecMtx, CondList] = ParseCondSpec(fid,nvariates,nconditions,nthcondspec)

  CondSpecMtx = [];
  CondList = [];

  ncondlist = readkeyvalue(fid,'nconditionlist',1,'%d');
  if(isempty(ncondlist)) return; end

  CondList = readkeyvalue(fid,'conditionlist',ncondlist,'%d');
  if(isempty(CondList)) return; end

  nsubblocks = readkeyvalue(fid,'nsubblocks',1,'%d');
  if(isempty(nsubblocks)) return; end

  for nthsubblock = 1:nsubblocks
    fgetnextline(fid);
    SubBlock = ParseSubBlock(fid,nvariates);
    if(isempty(SubBlock))
       fprintf('ERROR reading subblock %d from condspec %d in hs file %s\n',...
            nthsubblock,nthcondspec,fopen(fid));
       return;
    end
    CondSpecMtx = [CondSpecMtx SubBlock];

  end 

return

%-----------------------------------------------------------------------%
function SubBlock = ParseSubBlock(fid,nvariates)
  SubBlock = [];
    
  weight = readkeyvalue(fid,'weight',1,'%f');
  if(isempty(weight)) return; end

  subblocklen = readkeyvalue(fid,'length',1,'%f');
  if(isempty(subblocklen)) return; end

  method = readkeyvalue(fid,'method',1,'%s');
  if(isempty(method)) return; end

  switch(method)

    case 'replication'
      vect =  readkeyvalue(fid,'vector',subblocklen,'%f');
      vect = reshape(vect, [1 length(vect)]);
      if(isempty(vect)) return; end
      SubBlock = repmat(vect, [nvariates 1]);

    case 'shift'
      vect =  readkeyvalue(fid,'vector',subblocklen,'%f');
      vect = reshape(vect, [1 length(vect)]);
      if(isempty(vect)) return; end
      SubBlock = [];
      for n = 1:nvariates,
        v2 = fast_mshift(vect,[0 n-1],0);
	SubBlock = [SubBlock; v2];
      end

    case 'parametric'
      functionname = readkeyvalue(fid,'function',1,'%s');
      if(isempty(functionname)) return; end
      nparams = nparamsfromname(functionname);
      if(nparams < 0) return; end
      for n = 1:nvariates,
        parlist = readkeyvalue(fid,'parameters',nparams,'%f');
        if(isempty(parlist)) return; end
        vect = parfunction(functionname,parlist,subblocklen);
	SubBlock = [SubBlock; vect];
      end

    otherwise
      fprintf('ParseSubBlock: method %s unrecognized\n',method);
      return;

  end % switch(method)

  SubBlock = weight*SubBlock;

return;

%------------------------------------------------------------%
function nparams = nparamsfromname(parfunction)

  nparams = -1;

  switch(parfunction)
  
    case 'gamma', 
      nparams = 3;
      return;

    otherwise
      fprintf('Parameter function %s unrecognized\n',parfunction);
      return

  end % switch(parfunction)

return;
%------------------------------------------------------------%
function vect = parfunction(functionname,parlist,subblocklen)

  vect = [];

  switch(functionname)
  
    case 'gamma', 
      vect = zeros(subblocklen,1);
      delta = parlist(1);
      tau   = parlist(2);
      dt    = parlist(3);
      t = dt*[0:subblocklen-1];
      itlz = find(t<=0);
      itgz = find(t>0);
      vect(itlz) = 0;
      vect(itgz) = ((t(itgz)-delta)/tau).^2 .* exp(-(t(itgz)-delta)/tau);
      vect = reshape(vect, [1 subblocklen]);
      return;

    otherwise
      fprintf('Parameter function %s unrecognized\n',functionname);
      return

  end % switch(parfunction)

return

%------------------------------------------------------------%
function val = readkeyvalue(fid,keyexpected,nexpected,valfmt)

  val = [];
  line = fgetnextline(fid);
  if(isempty(line)) 
    fprintf('Error trying to read key %s\n',keyexpected);
    fprintf('  End-of-file reached\n');
    return; 
  end

  [key count errmsg nextind] = sscanf(line,'%s',1);
  if(~isempty(errmsg))
    fprintf('Error trying to read key %s\n',keyexpected);
    fprintf('%s\n',errmsg);
    return;
  end
  if(~strcmp(keyexpected,key))
    fprintf('Error trying to read key %s\n',keyexpected);
    fprintf('  actual key is %s\n',key);
    return;
  end

  % strip off key text %
  line = line(nextind:length(line));
  [val count errmsg ] = sscanf(line,valfmt,nexpected);
  if(~isempty(errmsg))
    fprintf('Error trying to read key %s\n',keyexpected);
    fprintf('%s\n',errmsg);
    return;
  end

return

%-------------------------------------------------------------------%
% Gets the next line that is not empty and does not begin with '#'.
% Converts the line to lower case
function line = fgetnextline(fid)
  line = [];
  while(~feof(fid) & (length(line) == 0 | line(1) == '#'))
    line = fgetl(fid);
  end
  lower(line);
return;

