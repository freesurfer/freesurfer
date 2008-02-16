function str = remove_spaces(str)
% str = remove_spaces(str)

for j=length(str):-1:1
		if (isspace(str(j)))
			 str(j) = 0 ;
		else
			 break ;
		end
end

