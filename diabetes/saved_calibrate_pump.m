
% A is active insulin
% M is measured blood glucose level
% F is carbs (in grams)
% D is insulin dose given
% I  - BG sensitivity to carbs (amount BG is raised/unit carb)



A =  [1.2  .4; .225+.65 .35; .225+.775 .35; .25+.725 .45];
F = [21 29 21 21];

M =  [202  136; 153 162; 119 152; 174 235; 190 255];
A =  [1.2  .4; .225+.65 .35; .65 .15; .225+.775 .35; .25+.725 .45];
F = [21 29 32 21 21];

%M = M(1:2,:);
%A = A(1:2,:);
%F = F(1:2);
dt = 130/60;  % in hours


if 0
M = [136 123 147 102 84];
F = [99 15 0 35 0];
dt = [0 120 113 18 58]/60;
A = [.4+2.35 .94 .275 .15 .8];
end


if 1 
% 3/7
M =  [76 106 89] ;
F =  [15 20 40] ;
A = [.6 .225 .3] ;
D = [0 .45 .9] ;
dt = [0  85 67] / 60;
% 3/6
M =  [201 228 235 279 204];
F =  [0   23   0   0   0]
A = [.6 .3 .6   .55 .925];  % active before dose
D = [0 1.1  .175  .8 0];     % dose given
dt = [0  87 78 30 39]/60;


end

S = 150 ;
e = .0 ;
I = 7;
C = 25;

min_rms = 10000;

for S=50:300
  disp(sprintf('processing S = %d', S)) 
  drawnow;
  for e=-.1:.025:.1
     for I=3:20
       for C=35:35
%       	   rms = compute_bg_error_pairs(M,A,F,dt,S,e,I,C);
       	   rms = compute_bg_error(M,A,D,F,dt,S,e,I,C);
	   if (rms < min_rms)
	      min_rms = rms ;
	      Sbest = S ;
	      ebest = e ;
	      Ibest = I ;
	      Cbest = C ;
	      disp(sprintf('new optimum %2.0f found at S = %d, e=%2.2f, I=%d, C=%d', rms,S,e,I,C)) 
	   end
	end
     end
  end
end

disp(sprintf('global optimum %2.0f found at S = %d, e=%2.2f, I=%d, C=%d', min_rms,Sbest,ebest,Ibest,Cbest)) ;

error = 0 ;

if 0
dt = 130/60;  % in hours
npairs = size(M,1);
for i=1:npairs
    insulin = A(i,1)-A(i,2) ;
    M1 = predict_bg(M(i,1), Sbest, insulin, ebest, dt, F(i), Cbest, Ibest);
    error = abs(M(i,2)-M1) ;
    disp(sprintf('tp %d, BG %d-->%d (%2.0f predicted), insulin=%2.2f, carbs=%d',i, M(i,1), M(i,2), M1,insulin,F(i)));
end

else
ntps = size(M,2);
for i=2:ntps
    insulin = A(i-1) + D(i-1)- A(i) ;
    M1 = predict_bg(M(i-1), Sbest, insulin, ebest, dt(i), F(i-1), Cbest, Ibest);
    error = abs(M(i)-M1) ;
    disp(sprintf('tp %d, BG %d-->%d (%2.0f predicted), insulin=%2.2f, carbs=%d, error=%2.0f',...
    		     i, M(i-1), M(i), M1,insulin,F(i-1),error));
end
end
