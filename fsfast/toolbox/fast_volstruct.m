function fv = fast_mrivolstruct
% fastvol = fast_mrivolstruct

fv.data      = [];
fv.szvol     = []; % not necesarily same as size(fv.data)
fv.szpix     = [];
fv.nvox      = 0;
fv.rows      = [];
fv.cols      = [];
fv.slices    = [];
fv.planes    = [];
fv.volid     = '';
fv.precision = '';
fv.endian    = '';
fv.format    = '';
