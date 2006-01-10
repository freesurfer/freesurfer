function [havg, hstd, Nh] = hsa_convert(hsa,Nnnc)
  Nc = Nnnc+1;

  if(size(hsa,2) == 1 & length(size(hsa)) == 2)
    ud.nRows = 1;
    ud.nCols = 1;
    Nch2 = size(hsa,1);
  else
    [ud.nRows ud.nCols Nch2] = size(hsa);
  end
  Nh = Nch2/(2*Nc);

  hsa2 = permute(hsa, [3 2 1 4]); % t,c,r
  hsa3 = reshape(hsa2, [Nh 2 Nc ud.nCols ud.nRows ]); % h stat cond col row
  clear hsa2;

  hsa4 = permute(hsa3, [5 4 3  1 2 ]); % row col cond h stat
  clear hsa3;

  havg = (hsa4(:,:,:,:,1)); % col row cond havg 
  hstd = (hsa4(:,:,:,:,2)); % col row cond hstd

  %havg = squeeze(hsa4(:,:,:,:,1)); % col row cond havg 
  %hstd = squeeze(hsa4(:,:,:,:,2)); % col row cond hstd
  clear hsa4;

return

