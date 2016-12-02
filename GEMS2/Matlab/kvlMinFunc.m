function [x,f,exitflag,output] = kvlMinFunc( funObj, x0 )
% kvlMinFunc( funObj, x0 )
%
% Stripped-down minFunc with the following options hard-coded:
%
%  options = [];
%  options.display = 'full';
%  options.Method = 'cg';
%  options.LS = 0;
%



verbose = 1;
verboseI= 1;
debug = 1;
doPlot = 0;

%
LS = 0;  % Type of line search
LS_init = 2;  % Type of initial step size guess
cgUpdate = 2; % Type of non-linear CG algorithm
%  cgUpdate = 0; % Type of non-linear CG algorithm - seems to go smoother down than 2
maxFunEvals = 1000;
maxIter = 500;
tolFun = 1e-5;
tolX = 1e-9;
c1 = 1e-4;


% Initialize
p = length(x0);
d = zeros(p,1);
x = x0;
t = 1;

% If necessary, form numerical differentiation functions
funEvalMultiplier = 1;


% Evaluate Initial Point
[f,g] = funObj(x);
funEvals = 1;



% Output Log
if verboseI
    fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
end



% Initialize Trace
trace.fval = f;
trace.funcCount = funEvals;

% Check optimality of initial point
if sum(abs(g)) <= tolFun
    exitflag=1;
    msg = 'Optimality Condition below TolFun';
    if verbose
        fprintf('%s\n',msg);
    end
    if nargout > 3
        output = struct('iterations',0,'funcCount',1,...
            'firstorderopt',sum(abs(g)),'message',msg,'trace',trace);
    end
    return;
end

% Perform up to a maximum of 'maxIter' descent steps:
for i = 1:maxIter

    % ****************** COMPUTE DESCENT DIRECTION *****************
    if i == 1
        d = -g; % Initially use steepest descent direction
    else
        gtgo = g'*g_old;
        gotgo = g_old'*g_old;

        if cgUpdate == 0
            % Fletcher-Reeves
            beta = (g'*g)/(gotgo);
        elseif cgUpdate == 1
            % Polak-Ribiere
            beta = (g'*(g-g_old)) /(gotgo);
        elseif cgUpdate == 2
            % Hestenes-Stiefel
            beta = (g'*(g-g_old))/((g-g_old)'*d);
        else
            % Gilbert-Nocedal
            beta_FR = (g'*(g-g_old)) /(gotgo);
            beta_PR = (g'*g-gtgo)/(gotgo);
            beta = max(-beta_FR,min(beta_PR,beta_FR));
        end

        d = -g + beta*d;

        % Restart if not a direction of sufficient descent
        if g'*d > -tolX
            if debug
                fprintf('Restarting CG\n');
                % figure
                % plot( trace.fval )
                % pause
            end
            beta = 0;
            d = -g;
        end

        % Old restart rule:
        %if beta < 0 || abs(gtgo)/(gotgo) >= 0.1 || g'*d >= 0

    end
    g_old = g;


    %  if ~isLegal(d)
    %      fprintf('Step direction is illegal!\n');
    %      pause;
    %      return
    %  end

    % ****************** COMPUTE STEP LENGTH ************************

    % Directional Derivative
    gtd = g'*d;

    % Check that progress can be made along direction
    if gtd > -tolX
        exitflag=2;
        msg = 'Directional Derivative below TolX';
        break;
    end

    % Select Initial Guess
    if i == 1
        t = min(1,1/sum(abs(g)));
    else
        if LS_init == 0
            % Newton step
            t = 1;
        elseif LS_init == 1
            % Close to previous step length
            t = t*min(2,(gtd_old)/(gtd));
        elseif LS_init == 2
            % Quadratic Initialization based on {f,g} and previous f
            t = min(1,2*(f-f_old)/(gtd));
        elseif LS_init == 3
            % Double previous step length
            t = min(1,t*2);
        end

        if t <= 0
            t = 1;
        end
    end
    f_old = f;
    gtd_old = gtd;

    fr = f;

    % Line Search
    f_old = f;
    % Use Armijo Bactracking
    % Perform Backtracking line search
    if 0
      [t,x,f,g,LSfunEvals] = ArmijoBacktrack(x,t,d,f,fr,g,gtd,c1,LS,tolX,debug,doPlot,1,funObj);
    elseif 1
      [t,x,f,g,LSfunEvals] = kvlBackTrack(x,t,d,f,fr,g,gtd,c1,LS,tolX,debug,doPlot,1,funObj);
    else
      % Somehow this should work the same as kvlBackTrack, but it doesn't ???


      % Evaluate the Objective and Gradient at the Initial Step
      [f_new,g_new] = funObj(x + t*d);
      LSfunEvals = 1;

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
          LSfunEvals = LSfunEvals+1;

          % Check whether step size has become too small
          if sum(abs(t*d)) <= tolX
              if debug
                  fprintf('Backtracking Line Search Failed\n');
              end
              t = 0;
              f_new = f;
              g_new = g;
              break;
          end
      end

      x_new = x + t*d;

      x = x_new;
      f = f_new;
      g = g_new;

      if ( LSfunEvals > 1 )
        disp( [ '          needed to do ' num2str( LSfunEvals-1 ) ' step halving trials' ] )
      end
    end

    funEvals = funEvals + LSfunEvals;


    % Output iteration information
    if verboseI
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,t,f,sum(abs(g)));
    end


    % Update Trace
    trace.fval(end+1,1) = f;
    trace.funcCount(end+1,1) = funEvals;

    % Check Optimality Condition
    if sum(abs(g)) <= tolFun
        exitflag=1;
        msg = 'Optimality Condition below TolFun';
        break;
    end

    % ******************* Check for lack of progress *******************

    if sum(abs(t*d)) <= tolX
        exitflag=2;
        msg = 'Step Size below TolX';
        break;
    end


    if abs(f-f_old) < tolX
        exitflag=2;
        msg = 'Function Value changing by less than TolX';
        break;
    end

    % ******** Check for going over iteration/evaluation limit *******************

    if funEvals*funEvalMultiplier > maxFunEvals
        exitflag = 0;
        msg = 'Exceeded Maximum Number of Function Evaluations';
        break;
    end

    if i == maxIter
        exitflag = 0;
        msg='Exceeded Maximum Number of Iterations';
        break;
    end

end

if verbose
    fprintf('%s\n',msg);
end
if nargout > 3
    output = struct('iterations',i,'funcCount',funEvals*funEvalMultiplier,...
        'firstorderopt',sum(abs(g)),'message',msg,'trace',trace);
end


