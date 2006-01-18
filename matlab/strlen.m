function len = strlen(str)
%  len = strlen(str)
% compute the # of characters in str (ignoring 0s at the end)

len = length(str) ;
for i=length(str):-1:1
	if (str(i) ~= 0)
		break ;
	end
	
	len = len-1;
end

		
