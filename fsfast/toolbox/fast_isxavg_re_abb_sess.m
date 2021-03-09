% fast_isxavg_re_abb_sess.m
%


%
% fast_isxavg_re_abb_sess.m
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

% These variable must be set first
%SessList = splitstring('$SessList');
%fsd      = '$fsd';
%analysis = '$analysis';
%hemi = splitstring('$hemi');
%spacedir = '$spacedir';
%outdir = '$outdir';
%usepct = [$usepct];
%synth = [$synth];

tic;

if(synth)
  st = sum(100*clock);
  randn('state',st);
  fprintf('Synthsizing %g\n',st);
end

nsess = size(SessList,1);
nhemi = size(hemi,1);

for nthhemi = 1:nhemi

  hid = deblank(hemi(nthhemi,:));
  if(strcmp(hid,'nohemi'))  
    hemicode = '';
  else                       
    fprintf('hemi = %s   (%g)\n',hid,toc);
    hemicode = sprintf('-%s',hid);
  end

  % get the nslices from the first sess
  sessdir = deblank(SessList(1,:));
  hstem = sprintf('%s/%s/%s/%s/h%s',sessdir,fsd,analysis,spacedir,hemicode);
  [nrows ncols nframes fs nslices endian bext] = fmri_bfiledim(hstem);
  if(isempty(nrows))
    fprintf('ERROR: loading %s\n',hstem);
    return;
  end

  % Loop over each slice separately %
  fprintf('slice ');
  for slice = 0:nslices-1
    fprintf('%d ',slice);

    % Loop over each session, collect data %
    %fprintf('Loading data\n');
    hrall = zeros(nsess,nrows*ncols);
    hiall = zeros(nsess,nrows*ncols);
    for nthsess = 1:nsess
      sessdir = deblank(SessList(nthsess,:));
      hstem = sprintf('%s/%s/%s/%s/h%s',sessdir,fsd,...
		      analysis,spacedir,hemicode);
      if(~synth)
        [h mristruct] = fast_ldbslice(hstem,slice);
        if(isempty(h))
          fprintf('ERROR: loading %s\n',hstem);
          return;
        end
	if(isempty(mristruct)) mristruct = fast_mri_struct; end
	mristruct.voldim = [ncols nrows nslices];
        hr = h(:,:,8);  % real
        hi = h(:,:,9);  % imag
      else
        mristruct = fast_ldbhdr(hstem);
        hr = randn(nrows,ncols);
        hi = randn(nrows,ncols);
      end
      if(usepct)
        h0 = h(:,:,11); % mean offset (for pct)
        ind = find(h0==0);
        h0(ind) = 10e10;
        hr = hr./h0;
        hi = hi./h0;
      end
      hrall(nthsess,:) = reshape1d(hr)';
      hiall(nthsess,:) = reshape1d(hi)';
    end % sess

    %fprintf('Analyzing\n');
    % Analyze the data %

    y = [hrall; hiall];
    X = zeros(2*nsess,2);
    X(1:nsess,1) = 1;
    X(nsess+1:end,2) = 1;

    %y = reshape(avball,[nrows*ncols nsess])';
    %X = ones(nsess,1);
    [beta, rvar, vdof] = fast_glmfit(y,X);

    for c = 1:3
     if(c==1)
       % cosine only
       C = [1 0]; 
       sigstem = sprintf('%s/cos/sig%s',outdir,hemicode);
     elseif(c==2)
       % sin only
       C = [0 1]; 
       sigstem = sprintf('%s/sin/sig%s',outdir,hemicode);
     else
       C = eye(2); % avb
       sigstem = sprintf('%s/avb/sig%s',outdir,hemicode);
     end

     [F, Fsig, ces] = fast_fratio(beta,X,rvar,C);
     if(c == 3)
       %[m i] = max(abs(ces));
       %ind = sub2ind(size(ces),i,1:nrows*ncols);
       sig = -sign(ces(2,:)).*log10(abs(Fsig)); % give sign of sin
     else
       sig = -sign(ces).*log10(abs(Fsig));
     end
     sig = reshape(sig,[nrows ncols]);
     fast_svbslice(sig,sigstem,slice,'bfloat',mristruct);
    end

  end % slice
  fprintf('\n');

end % hemi

fprintf('done %g\n',toc);

