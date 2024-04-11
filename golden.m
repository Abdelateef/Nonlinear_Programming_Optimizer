% Golden Section Optimization Method
% Inputs:
%   xl: Lower bound of the initial interval
%   xu: Upper bound of the initial interval
%   error: Tolerance for the convergence criterion
%   f: Objective function to be minimized or maximized
%   option: 'Min' for minimization, 'Max' for maximization
function golden(xl,xu,error,f,option)
n=log(error/(xu-xl))/log(0.618);
n=ceil(n);
c=0;
for i=1:n
    d=0.618*(xu-xl);
    x1=xl+d;
    x2=xu-d;
    if (option=='Max')
        if (f(x1)<=f(x2))
            xu=x1;
        else
            xl=x2;
        end 
    end
    if (option=='Min')
        if (f(x1)>=f(x2))
            xu=x1;
        else
            xl=x2;
        end 
    end 
    c=c+1;
end
fprintf('--------------------The Golden Section Method--------------------\n');
fprintf('The final interval at which function is %s [%f , %f]\n',option,xl,xu);
fprintf('the %s value is %f and it functon value is %f\n',option,((xl+xu)/2),f(((xl+xu)/2)));
end 
