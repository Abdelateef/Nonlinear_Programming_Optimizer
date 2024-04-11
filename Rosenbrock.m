% Rosenbrock Function along a Direction
% Inputs:
%   X1, X2: Initial values for the Rosenbrock function
% Output:
%   ff_lamda: Function representing the Rosenbrock function along a direction

function ff_lamda = Rosenbrock(X1,X2)
syms x_1 x_2 lamda;
f= @(x1,x2) 100*(x2-x1^2)^2+(1-x1)^2;
f_sym=100*(x_2-x_1^2)^2+(1-x_1)^2;
gradient_sym = gradient(f_sym, [x_1, x_2]);
s = subs(gradient_sym, [x_1, x_2], [0, 0]);
lamdas_sub= [X1;X2]+lamda*s;
ff= f(lamdas_sub(1),lamdas_sub(2));
ff_lamda = matlabFunction(ff, 'Vars', lamda);
end
