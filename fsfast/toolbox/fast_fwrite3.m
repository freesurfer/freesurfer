function fast_fwrite3(fid, val)
% fast_fwrite3(fid, val)
%
% Write a 3 byte integer into a file
%

if(nargin ~= 2)
  fprintf('USAGE: fast_fwrite3(fid, val)\n');
  return;
end

%fwrite(fid, val, '3*uchar') ;
b1 = bitand(bitshift(val, -16), 255) ;
b2 = bitand(bitshift(val, -8), 255) ;
b3 = bitand(val, 255) ; 
fwrite(fid, b1, 'uchar') ;
fwrite(fid, b2, 'uchar') ;
fwrite(fid, b3, 'uchar') ;

return
