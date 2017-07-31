

for delta = di:di:insulin_init(2)
    insulin_tmp = insulin_init;
    insulin_tmp(2) = insulin_init(2) - di ;
    insulin_tmp(3) = insulin_init(2) + di ;
end

