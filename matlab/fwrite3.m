function fwrite3(fid, val)

% write a 3 byte integer out of a file
%fwrite(fid, val, '3*uchar') ;
b1 = bitand(bitshift(val, -16), 255) ;
b2 = bitand(bitshift(val, -8), 255) ;
b3 = bitand(val, 255) ; 
fwrite(fid, b1, 'uchar') ;
fwrite(fid, b2, 'uchar') ;
fwrite(fid, b3, 'uchar') ;


