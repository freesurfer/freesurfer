function fmri_svcvm(cvmstruct,stem)
% fmri_svcvm(cvmstruct,stem)

if(nargin ~= 2)
     msg = 'USAGE: fmri_svcvm(cvmstruct,stem)';
     qoe(msg); error(msg);
end

cvmstruct.sz = size(cvmstruct.cvm,1);

fname = sprintf('%s.bfloat',stem);
fmri_svbfile(cvmstruct.cvm,fname);
fname = sprintf('%s.cvm',stem);
fid = fopen(fname,'w');
if(fid == -1) 
  msg = sprintf('Could not open %s for writing',fname);
  qoe(msg);error(msg);
end
fprintf(fid,'%d\n',cvmstruct.version);
fprintf(fid,'%d\n',cvmstruct.n);
fprintf(fid,'%f\n',cvmstruct.d);
fprintf(fid,'%d\n',cvmstruct.sz);
fprintf(fid,'%d\n',cvmstruct.norm);
fprintf(fid,'%d\n',cvmstruct.inv);
if(~isempty(cvmstruct.acoravg) & ~isempty(cvmstruct.acorstd))  
  fprintf(fid,'%f\n',cvmstruct.acoravg);
  fprintf(fid,'%f\n',cvmstruct.acorstd);
end

fclose(fid);

return;
