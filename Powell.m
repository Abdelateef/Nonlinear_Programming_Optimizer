% Powell Function along a Direction
% Inputs:
%   X1, X2, X3, X4: Initial values for the Powell function
% Output:
%   ff_lamda: Function representing the Powell function along a direction

function ff_lamda = Powell(X1,X2,X3,X4)
syms x_1 x_2 x_3 x_4 lamda;
f= @(x1,x2,x3,x4) (x1+10*x2)^2 + 5*(x3-x4)^2 +(x2-2*x3)^4 + 10*(x1-x4)^4;
f_sym= (x_1+10*x_2)^2 + 5*(x_3-x_4)^2 +(x_2-2*x_3)^4 + 10*(x_1-x_4)^4;
gradient_sym = gradient(f_sym, [x_1, x_2,x_3,x_4]);
s = subs(gradient_sym, [x_1, x_2,x_3,x_4], [1,0,0,0]);
lamdas_sub= [X1;X2;X3;X4]+lamda*s;
ff= f(lamdas_sub(1),lamdas_sub(2),lamdas_sub(3),lamdas_sub(4));
ff_lamda = matlabFunction(ff, 'Vars', lamda);
end
