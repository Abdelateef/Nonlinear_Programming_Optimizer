% Fibonacci Optimization Method
% Inputs:
%   a: Lower bound of the initial interval
%   b: Upper bound of the initial interval
%   error: Tolerance for the convergence criterion
%   f: Objective function to be minimized or maximized
%   option: 'Min' for minimization, 'Max' for maximization
function fib(a,b,error,f,option)
fnn=(b-a)/error;
fibonacciArray = [0, 1];
while true
    nextValue = fibonacciArray(end) + fibonacciArray(end-1);
    fibonacciArray = [fibonacciArray, nextValue];
    if (fibonacciArray(end) + fibonacciArray(end-1) >= fnn)
        nextValue = fibonacciArray(end) + fibonacciArray(end-1);
        fibonacciArray = [fibonacciArray, nextValue];
        rk=[]; 
        ennd=length(fibonacciArray); 
        for i=3:ennd
            rk=[rk,(fibonacciArray(ennd-1)/fibonacciArray(ennd))];
            ennd=ennd-1;
        end
        break;
    end
end
% the goal to get the min
c=0;
for i=1:length(rk)
    x1=a+rk(i)*(b-a);
    x2=b-rk(i)*(b-a);
    if (option=='Max')
        if (f(x1)<=f(x2))
            b=x1;
        else
            a=x2;
        end 
    end
    if (option=='Min')
        if (f(x1)>=f(x2))
            b=x1;
        else
            a=x2;
        end 
    end 
    c=c+1;
end
fprintf('--------------------The Fibonacci Method------------------------\n');
fprintf('The final interval at which function is %s [%f , %f]\n',option,a,b);
fprintf('the %s value is %f and it functon value is %f\n',option,((a+b)/2),f(((a+b)/2)));
end  