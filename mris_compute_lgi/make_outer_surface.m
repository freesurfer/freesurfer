function make_outer_surface (filled_volume, se_diameter, output_surface)

% Original Author: Marie Schaer
% Date: 2007/11/14
% 
% This function takes as an input the binary volume resulting of the
% filling of the any surface (usually the pial one) using mris_fill, and
% will close the sulci using morphological operation, with a sphere as the
% structural element.
% 
% Parameters: 
% se_diameter is the diameter of the sphere (in mm), use 15mm by default 
% to close the sulci.  
% 
% Utilities to write the surface to the freesurfer format are modified from
% "freesurfer_write_surf" (from Darren Weber's bioelectromagnetism toolbox), 
% according to a suggestion by Don Hagler in FreeSurfer's mailing list on 
% August 3, 2007.
%
% Example: make_outer_surface('lh.pial.mgz',15,'lh.outer-pial')


    fprintf('reading filled volume...\n');
    vol=MRIread(filled_volume);
    volume=vol.vol;
    
    fprintf('closing volume...\n');
    se=strel('ball',se_diameter,se_diameter);
    BW2=imclose(volume,se);
    thresh = max(BW2(:))/2;
    i=find(BW2<=thresh);
    BW2(i)=0;
    i=find(BW2>thresh);
    BW2(i)=255;
    
    [f,v] = isosurface(BW2,100);

    v2=[129-v(:,1) v(:,3)-129 129-v(:,2)]; % in order to cope with the different orientation 
    v=v2;

    fprintf('morphological closing done.\n');
    fprintf('writing outer surface...\n');
    
    fname=output_surface;
    vert = v;
    face = f - 1;
    vnum = size(vert,1);  
    fnum = size(face,1);  
    
    % open it as a big-endian file
    fid = fopen(fname, 'wb', 'b');
    TRIANGLE_FILE_MAGIC_NUMBER = 16777214;
    fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);
    
    % Ouput a couple of text lines with creation date
    str = sprintf('created from matlab on %s\n',datestr(now))
    fwrite(fid, str,'char');
    fwrite(fid, vnum,'int32');
    fwrite(fid, fnum,'int32');

    % reshape vert into column array and write
    vert = reshape(vert',size(vert,1)*size(vert,2),1);
    fwrite(fid, vert,'float32');

    % reshape face into column array and write
    face = reshape(face',size(face,1)*size(face,2),1);
    fwrite(fid, face,'int32');
    fclose(fid) ;
return
end
