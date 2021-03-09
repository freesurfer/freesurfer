function err = write_analyze_hdr(hdr,hdrfile)
% err = write_analyze_hdr(hdr,hdrfile)
%


%
% write_analyze_hdr.m
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



err = 1;

if(nargin ~= 2)
  fprintf('err = write_analyze_hdr(hdr,hdrfile)\n');
  return;
end

% Opening as big or little does not seem to matter
fp = fopen(hdrfile,'wb');
if(fp == -1) 
  fprintf('ERROR: could not open %s for writing.\n',hdrfile);
  return;
end

% Should do a better job of checking lengths

fwrite(fp,hdr.key.sizeof_hdr(1),'int'); % should always be 348
fwrite(fp,hdr.key.data_type(1:10),'char');
fwrite(fp,hdr.key.db_name(1:18),'char');
fwrite(fp,hdr.key.extents(1),'int');
fwrite(fp,hdr.key.session_error(1),'short');
fwrite(fp,hdr.key.regular(1),'char');
fwrite(fp,hdr.key.hkey_un0(1),'char');

fwrite(fp,hdr.dime.dim(1:8),'short');
fwrite(fp,hdr.dime.vox_units(1:4),'char');
fwrite(fp,hdr.dime.cal_units(1:8),'char');
fwrite(fp,hdr.dime.unused1,'short');
fwrite(fp,hdr.dime.datatype,'short');
fwrite(fp,hdr.dime.bitpix,'short');
fwrite(fp,hdr.dime.dim_un0,'short');
fwrite(fp,hdr.dime.pixdim(1:8),'float');
fwrite(fp,hdr.dime.vox_offset,'float');
fwrite(fp,hdr.dime.roi_scale,'float');
fwrite(fp,hdr.dime.funused1,'float');
fwrite(fp,hdr.dime.funused2,'float');
fwrite(fp,hdr.dime.cal_max,'float');
fwrite(fp,hdr.dime.cal_min,'float');
fwrite(fp,hdr.dime.compressed,'int');
fwrite(fp,hdr.dime.verified,'int');
fwrite(fp,hdr.dime.glmax,'int');
fwrite(fp,hdr.dime.glmin,'int');

fwrite(fp,hdr.hist.descrip(1:80),'char');
fwrite(fp,hdr.hist.aux_file(1:24),'char');
fwrite(fp,hdr.hist.orient,'char');
fwrite(fp,hdr.hist.originator(1:10),'unsigned char');
fwrite(fp,hdr.hist.generated(1:10),'char');
fwrite(fp,hdr.hist.scannum(1:10),'char');
fwrite(fp,hdr.hist.patient_id(1:10),'char');
fwrite(fp,hdr.hist.exp_date(1:10),'char');
fwrite(fp,hdr.hist.exp_time(1:10),'char');
fwrite(fp,hdr.hist.hist_un0(1:3),'char');
fwrite(fp,hdr.hist.views,'int');
fwrite(fp,hdr.hist.vols_added,'int');
fwrite(fp,hdr.hist.start_field,'int');
fwrite(fp,hdr.hist.field_skip,'int');
fwrite(fp,hdr.hist.omax,'int');
fwrite(fp,hdr.hist.omin,'int');
fwrite(fp,hdr.hist.smax,'int');
fwrite(fp,hdr.hist.smin,'int');
fclose(fp);

% Save vox2ras (if there) as mat file 
if(~isempty(hdr.vox2ras))
  basename = hdrfile(1:end-4);
  matfile = sprintf('%s.mat',basename);
  M = hdr.vox2ras;
  save(matfile,'M','-v4');
end

return;





