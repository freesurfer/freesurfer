function [avglist, stdlist] = fast_hlist(Nc,Nh)
% [avglist stdlist] = fast_hlist(Nc,Nh)

avglist = [];
stdlist = [];

if(nargin ~= 2)
  msg = '[avglist stdlist] = fast_hlist(Nc,Nh)';
  qoe(msg); error(msg);
end 

n = 1;
for c = 1:Nc
  for s = 1:2
    for h = 1:Nh
      if(s==1) avglist = [avglist n];
      else     stdlist = [stdlist n];
      end
      n = n + 1;
    end
  end
end

return;
