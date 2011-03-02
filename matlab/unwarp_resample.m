function [imvol_out, M_out] = unwarp_resample(imvol,M,imvol_out_size,M_out,Mdc,unwarpflag,Jacobianflag,plotflag,interp_method,inflag,thruflag,gradfilename)

%
% In this file:
%
% unwarp_resample - from AD's unwarp_and_resample_vol
% proj, jdproj - for projecting out e.g. throughplane component of
%                displacement and jacobian correction
% interp_trilin, rcs2index - for interpolating in the
%                            gradwarpvolumes


%
% unwarp_resample.m
%
% Original Author: Elizabeth Haley
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.7 $
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



devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir); % for qoe()
if(isempty(findstr(d,path))) path(path,d); end

[nrows,ncols,nslices] = size(imvol);
nrows_out = imvol_out_size(1);
ncols_out = imvol_out_size(2);
nslices_out = imvol_out_size(3);
Mi = inv(M);
imvol_out = zeros(imvol_out_size);

if inflag & thruflag
  qoe('unwarp_and_resample_vol() should not be invoked with inflag and thruflag.');
  error('unwarp_and_resample_vol() should not be invoked with inflag and thruflag.');
end

if unwarpflag %E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if inflag
    fprintf(1,'Inplane-only unwarping...\n');
    unwarpdims=2;
    SGflag=0;
  elseif thruflag
    fprintf(1,'Throughplane-only unwarping...\n');
    unwarpdims=1;
    SGflag=0;
  else
    fprintf(1,'Full 3D unwarping...\n');
    unwarpdims=3;
    SGflag=1;
  end

  % setting SGflag based on unwarpdims is kind of a kludge - can't
  % really assume you'd never unwarp GE 3d.  Really, should check
  % gradfilename against some list you'd make, to confirm format.

  disp(sprintf('gradfilename is %s.  Loading...', gradfilename));

  fid = fopen(gradfilename,'rb','b');
  if(fid == -1)
    msg = sprintf('could not open %s\n',gradfilename);
    qoe(msg);
    error(msg);
  end
  ncoords_x = fread(fid,1,'int32');
  ncoords_y = fread(fid,1,'int32');
  ncoords_z = fread(fid,1,'int32');
  nv = ncoords_x*ncoords_y*ncoords_z;
  coords_x = fread(fid,2,'float32')*1000;
  coords_y = fread(fid,2,'float32')*1000;
  coords_z = fread(fid,2,'float32')*1000;
  Dx = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')'*1000;
  if(length(Dx) ~= nv)
    msg = sprintf('ERROR: reading Dx from %s, read %d, expected %d\n',...
		  gradfilename,nv,length(Dx));
    qoe(msg); error(msg);
  end
  Dy = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')'*1000;
  if(length(Dy) ~= nv)
    msg = sprintf('ERROR: reading Dy from %s, read %d, expected %d\n',...
		  gradfilename,nv,length(Dy));
    qoe(msg);  error(msg);
  end
  Dz = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')'*1000;
  if(length(Dz) ~= nv)
    msg = sprintf('ERROR: reading Dz from %s, read %d, expected %d\n',...
		  gradfilename,nv,length(Dz));
    qoe(msg);  error(msg);
  end


  if SGflag % Siemens
    JacDet = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')';
    % Should do a check here
  else %GE
    JacDetx = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')';
    JacDety = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')';
    JacDetz = fread(fid,ncoords_x*ncoords_y*ncoords_z,'float32')';
  end
    
  fclose(fid);
  fprintf(1,'Unwarp data loaded.\n');

  % maxjd = 3; % moved Jacobian density correction cutoff to jdproj()
  Mgw = [0,(coords_x(2)-coords_x(1))/(ncoords_x-1),0,coords_x(1)-(coords_x(2)-coords_x(1))/(ncoords_x-1);
	  (coords_y(2)-coords_y(1))/(ncoords_y-1),0,0,coords_y(1)-(coords_y(2)-coords_y(1))/(ncoords_y-1);
	  0,0,(coords_z(2)-coords_z(1))/(ncoords_z-1),coords_z(1)-(coords_z(2)-coords_z(1))/(ncoords_z-1);
	  0, 0, 0, 1];
  Mgwi = inv(Mgw); 
end %if unwarpflag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if(plotflag) figure(1); end

tic;
for s = [1:nslices_out]
%  if (mod(s,10)==0)
    fprintf('slice %3d/%d  %g\t\n',s,nslices_out,toc)
%  end

  % Get CRS in target volume %
  [R,C] = meshgrid(1:nrows_out,1:ncols_out);
  R = reshape(R,1,nrows_out*ncols_out);
  C = reshape(C,1,nrows_out*ncols_out);
  S = ones(size(C))*s;
  % Indices corresponding to the CR subscripts in the slice %
  imindx = R+(C-1)*nrows_out;

  % Convert CRS to XYZ in target volume %
  X = M_out(1,1)*R + M_out(1,2)*C + M_out(1,3)*S + M_out(1,4);
  Y = M_out(2,1)*R + M_out(2,2)*C + M_out(2,3)*S + M_out(2,4);
  Z = M_out(3,1)*R + M_out(3,2)*C + M_out(3,3)*S + M_out(3,4);

  if ~unwarpflag
    % Then this is the correct XYZ, proceed %
    X0 = X;
    Y0 = Y;
    Z0 = Z;
    
    % Convert XYZ to CRS in source volume %
    R0 = Mi(1,1)*X0 + Mi(1,2)*Y0 + Mi(1,3)*Z0 + Mi(1,4);
    C0 = Mi(2,1)*X0 + Mi(2,2)*Y0 + Mi(2,3)*Z0 + Mi(2,4);
    S0 = Mi(3,1)*X0 + Mi(3,2)*Y0 + Mi(3,3)*Z0 + Mi(3,4);
    
    % Get the value at these locations from the input volume %
    im0vec = interp3(imvol,C0,R0,S0,interp_method);
    im0vec(find(isnan(im0vec))) = 0;
    clear im0vec2;
    im0vec2(imindx) = im0vec;
    im0 = reshape(im0vec2,nrows_out,ncols_out);
    imvol_out(:,:,s) = im0;
  else
    % Convert to XYZ in target volume to CRS in warp table %
    Rgw = Mgwi(1,1)*X + Mgwi(1,2)*Y + Mgwi(1,3)*Z + Mgwi(1,4);
    Cgw = Mgwi(2,1)*X + Mgwi(2,2)*Y + Mgwi(2,3)*Z + Mgwi(2,4);
    Sgw = Mgwi(3,1)*X + Mgwi(3,2)*Y + Mgwi(3,3)*Z + Mgwi(3,4);

    % Get the offset and project%
    Dx1 = interp_trilin(Rgw,Cgw,Sgw,Dx,ncoords_y,ncoords_x,ncoords_z);
    Dy1 = interp_trilin(Rgw,Cgw,Sgw,Dy,ncoords_y,ncoords_x,ncoords_z);
    Dz1 = interp_trilin(Rgw,Cgw,Sgw,Dz,ncoords_y,ncoords_x,ncoords_z);
    [Dx2,Dy2,Dz2] = proj(unwarpdims,Mdc,[Dx1;Dy1;Dz1]);


    if Jacobianflag
      if SGflag %Siemens
	JacDet1 = interp_trilin(Rgw,Cgw,Sgw,JacDet,ncoords_y,ncoords_x,ncoords_z);
      else %GE
	JacDetx1 = interp_trilin(Rgw,Cgw,Sgw,JacDetx,ncoords_y,ncoords_x,ncoords_z);
	JacDety1 = interp_trilin(Rgw,Cgw,Sgw,JacDety,ncoords_y,ncoords_x,ncoords_z);
	JacDetz1 = interp_trilin(Rgw,Cgw,Sgw,JacDetz,ncoords_y,ncoords_x,ncoords_z);
	JacDet1 = jdproj([JacDetx1;JacDety1;JacDetz1],Mdc,unwarpdims,3);
	% that last arg to jdproj is maxjd=3
	% The first arg is 3 rows x ncoords^3 columns
	% and JacDet1 is a row vector with ncoords^3 rows.
	
	% a = max(JacDet1 - JacDetx1.*JacDety1.*JacDetz1); keyboard;
	% Test: 3D, a==0, good.
      end
    end

    % Recompute the XYZ %
    X1 = X+Dx2;
    Y1 = Y+Dy2;
    Z1 = Z+Dz2;

    % Compute the CRS in source volume %
    R1 = Mi(1,1)*X1 + Mi(1,2)*Y1 + Mi(1,3)*Z1 + Mi(1,4);
    C1 = Mi(2,1)*X1 + Mi(2,2)*Y1 + Mi(2,3)*Z1 + Mi(2,4);
    S1 = Mi(3,1)*X1 + Mi(3,2)*Y1 + Mi(3,3)*Z1 + Mi(3,4);
  
    im1vec = interp3(imvol,C1,R1,S1,interp_method);
    im1vec(find(isnan(im1vec))) = 0;

    % test - reorder these before multiplying, too - ?
    % But they work this order for DG so I must have compensating errors.
    % No they didn't work for DG.
    clear im1vec2; %first mention
    im1vec2(imindx) = im1vec;

    if Jacobianflag
      imjdvec(imindx) = JacDet1; % a reordering of col<->row
      % imjdvec(1) == JacDet1(1), imjdvec(257) == JacDet1(2), ...
      im1vec2 = im1vec2.*imjdvec; % correct for voxel size ratio
      imjd = reshape(imjdvec,nrows_out,ncols_out); % for plot
    end
    im1 = reshape(im1vec2,nrows_out,ncols_out);
    imvol_out(:,:,s) = im1;
  end % if ~unwarpflag
  
%%   if plotflag
%%     figure(s);
%%     subplot(2,2,1); imagesc(im0); axis equal; colorbar; title(sprintf('Slice %d:  Original',s));
%%     if unwarpflag
%%       subplot(2,2,2); imagesc(im1); axis equal; colorbar; title(sprintf('Slice %d:  Unwarped',s));
%%       if Jacobianflag
%%         subplot(2,2,3); imagesc(imjd); axis equal; colorbar; title(sprintf('Slice %d:  Jacobian factor',s));
%%       end
%%     end
%%   end

end % for s

imvol_out(find(imvol_out<0))=0; % Clamp to positive values

return % u_r


%-------------------------------------------------%
function [dx,dy,dz] = proj(unwarpdims,Mdc,Din)

% Din is 3 rows by ncoords^3 columns, all ready to left-multiply by
% these 3x3's.  Dx,Dy,Dz is returned - each a vector('?) with
% ncoords^3 rows

% Even if we have C/R switched, that shouldn't affect this.

Mdci = inv(Mdc);
if unwarpdims == 1 % thruplane
  Dout = Mdc * diag([0;0;1]) * Mdci * Din; % keep slice-direction piece
elseif unwarpdims == 2 % inplane
  Dout = Mdc * diag([1;1;0]) * Mdci * Din; % keep column and row piece
else % full unwarp 
  Dout = Din;
end
dx = Dout(1,:);
dy = Dout(2,:);
dz = Dout(3,:);
return % proj

%-------------------------------------------------%
function JDout = jdproj(JDin,Mdc,unwarpdims,maxjd)

% JDin is 3 rows by ncoords^3 columns, all ready to left-multiply by
% these 3x3's.  JDout is returned, the row vector of ncoords^3 scalar
% jac det factors

% Even if we have C/R switched, that shouldn't affect this.

Mdci = inv(Mdc);
if unwarpdims == 1 % thruplane
  % thruplane: keep slice piece =JDslice, but set row/col piece to 1
  JDtmp = Mdc * diag([0;0;1]) * Mdci * JDin + ...
      Mdc * diag([1;1;0]) * Mdci * ones(3,size(JDin,2));
elseif unwarpdims == 2 % inplane
  % inplane: keep column and row piece =JDrow.*JDcol; set slice piece to 1
  JDtmp = Mdc * diag([1;1;0]) * Mdci * JDin + ...
      Mdc * diag([0;0;1]) * Mdci * ones(3,size(JDin,2));
else % full unwarp - multiply all three of JDrow .* JDcol .* JDslice
  JDtmp = JDin;
end
% derive a scalar from each column - to return a 1-row ncoords^3-cols object 
JDout = JDtmp(1,:) .* JDtmp(2,:) .* JDtmp(3,:);
JDout(find(JDout<0)) = 0; % or 1/maxjd? In my experience they're in [0.6,1.3]
JDout(find(JDout>maxjd)) = maxjd;
return % jdproj

%-------------------------------------------------%

function V_out = interp_trilin(R,C,S,M,nr,nc,ns)

devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir); % for qoe
if(isempty(findstr(d,path))) path(path,d); end

if(nargin ~= 7)
  msg = 'USAGE: V_out = interp_trilin(R,C,S,M,nr,nc,ns)';
  qoe(msg);  error(msg);
end

R = max(1,min(nr-1,R));
C = max(1,min(nc-1,C));
S = max(1,min(ns-1,S));
Rint = floor(R);
Cint = floor(C);
Sint = floor(S);
Rfrac = R-Rint;
Cfrac = C-Cint;
Sfrac = S-Sint;
%V_out = M(rcs2index(round(R),round(C),round(S),nr,nc,ns));

Wr = [1-Rfrac;Rfrac];
Wc = [1-Cfrac;Cfrac];
Ws = [1-Sfrac;Sfrac];
V_out = 0;
for dr=[0,1]
  for dc=[0,1]
    for ds=[0,1]
      V_out = V_out+M(rcs2index(Rint+dr,Cint+dc,Sint+ds,nr,nc,ns)).*Wr(1+dr,:).*Wc(1+dc,:).*Ws(1+ds,:);
    end
  end
end

%V_out = M(rcs2index(Rint,Cint,Sint,nr,nc,ns));

return;


%----------------------------------------------------%
function indexvec = rcs2index(R,C,S,nr,nc,ns)

devdir = getenv('DEV');
d = sprintf('%s/fsfast/toolbox',devdir);
if(isempty(findstr(d,path))) path(path,d); end

if(nargin ~= 6)
  msg = 'USAGE: indexvec = rcs2index(R,C,S,nr,nc,ns)';
  qoe(msg);  error(msg);
end

indexvec = R+(C-1)*nr+(S-1)*nc*nr;
%indexvec = C+(R-1)*nc+(S-1)*nc*nr;

return;
%----------------------------------------------------%

