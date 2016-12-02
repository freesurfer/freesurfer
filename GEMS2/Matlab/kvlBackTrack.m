function  [t,x_new,f_new,g_new,funEvals,H] = kvlBackTrack(...
    x,t,d,f,fr,g,gtd,c1,LS,tolX,debug,doPlot,saveHessianComp,funObj,varargin)
%

% Evaluate the Objective and Gradient at the Initial Step
[f_new,g_new] = funObj(x + t*d);
funEvals = 1;

while f_new > fr + c1*t*gtd

    temp = t;
    % Backtrack w/ fixed backtracking rate
    %if debug
    %    fprintf('Fixed BT\n');
    %end
    t = 0.5*t;

    f_prev = f_new;
    t_prev = temp;
    [f_new,g_new] = funObj(x + t*d);
    funEvals = funEvals+1;

    % Check whether step size has become too small
    if sum(abs(t*d)) <= tolX
        if debug
            fprintf('Backtracking Line Search Failed\n');
            disp( '%%%%%%%%%%%%%%%%%%%%%%' )
            pause
        end
        t = 0;
        f_new = f;
        g_new = g;
        break;
    end
end

x_new = x + t*d;

if ( funEvals > 1 )
  disp( [ '          needed to do ' num2str( funEvals-1 ) ' step halving trials' ] )
end
