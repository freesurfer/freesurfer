function str = remove_spaces(str)
% str = remove_spaces(str)

for j=length(str):-1:1
		if (isspace(str(1,j)))
			 str(j) = [];
		else
			 break ;
		end
end

