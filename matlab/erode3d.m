function vout =  erode3d(vin, num_iter)
% vout = function erode3d(vin, num_iter)


[width, height, depth] = size(vin);
for n=1:num_iter
  vout = vin ;
  for x=2:width-1
    for y=2:height-1
      for z=2:depth-1
        vals = vin(x-1:x+1, y-1:y+1, z-1:z+1) ;
        vout(x,y,z) = min(vals(:));
      end
    end
  end
      vin = vout ;
end
