function dot = fmri_lddot(dotfile)
% dot = fmri_lddot(dotfile)

if(nargin ~= 1)
  msg = 'USAGE: dot = fmri_lddot(dotfile)';
  qoe(msg);error(msg);
end

dot = load(dotfile);
[n1 n2] = size(dot);
ntp = n1/2;
nsensors = n2;

dot = reshape(dot, [ntp 2 nsensors]);
dot = permute(dot, [2 3 1]);

return;

