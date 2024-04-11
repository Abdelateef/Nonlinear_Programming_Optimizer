% Fletcher-Reeves Conjugate Gradient Method for Optimization
% Inputs:
%   eR: Tolerance for Rosenbrock function
%   eP: Tolerance for Powell function
%   R_X0: Initial guess for Rosenbrock function
%   P_X0: Initial guess for Powell function
% Output:
%   None (Displays the final values at which functions are minimized)
function Fletcher_Reeves(eR,eP,R_X0,P_X0)
syms x_1 x_2 x_3 x_4 lamda;
%% Rosenbrock s parabolic using Fletcher
x__1=R_X0;
f_R = @(x1,x2) 100*(x2-x1^2)^2+(1-x1)^2;
f_R_sym=100*(x_2-x_1^2)^2+(1-x_1)^2;
R_gradient_sym = gradient(f_R_sym, [x_1, x_2]);
R_hessian_sym = hessian(f_R_sym, [x_1, x_2]);
%inital 
g_f1 = subs(R_gradient_sym, [x_1, x_2],[x__1(1),x__1(2)]);
A1 = subs(R_hessian_sym, [x_1, x_2], [x__1(1),x__1(2)]);
s1 = -subs(R_gradient_sym, [x_1, x_2],[x__1(1),x__1(2)]);
lamda_star1=(transpose(g_f1)*g_f1)/(transpose(s1)*A1*s1);
x__2=x__1+lamda_star1*s1;
g_f2 = subs(R_gradient_sym, [x_1, x_2],[x__2(1),x__2(2)]);
grad_norm=norm(g_f2, Inf);
k=2;
while (grad_norm>eR)
    A2 = subs(R_hessian_sym, [x_1, x_2], [x__2(1),x__2(2)]);
    beta2=(transpose(g_f2)*g_f2)/((transpose(g_f1)*g_f1));
    s2=-g_f2+beta2*s1;
    lamda_star2=(transpose(g_f2)*g_f2)/(transpose(s2)*A2*s2);
    x__3=x__2+lamda_star2*s2;
    x__2=x__3;
    g_f2 = subs(R_gradient_sym, [x_1, x_2],[x__2(1),x__2(2)]);
    grad_norm_prev=grad_norm;
    grad_norm=norm(g_f2, Inf);
    k=k+1;
    if (grad_norm>grad_norm_prev)
        break;
    end 
end 
fprintf('---------------The Rosenbrock’s parabolic using Reeves CG Method----------------\n');
fprintf('The final value at which function is min [%f,%f]\n',x__2(1),x__2(2));

%% Powell s quartic using Fletcher
x__1=P_X0;
f_P = @(x1,x2,x3,x4) (x1+10*x2)^2 + 5*(x3-x4)^2 +(x2-2*x3)^4 + 10*(x1-x4)^4;
f_P_sym=(x_1+10*x_2)^2 + 5*(x_3-x_4)^2 +(x_2-2*x_3)^4 + 10*(x_1-x_4)^4;
P_gradient_sym = gradient(f_P_sym, [x_1, x_2,x_3,x_4]);
P_hessian_sym = hessian(f_P_sym, [x_1, x_2,x_3,x_4]);

%inital 
g_f1 = subs(P_gradient_sym, [x_1,x_2,x_3,x_4],[x__1(1),x__1(2),x__1(3),x__1(4)]);
A1 = subs(P_hessian_sym, [x_1, x_2,x_3,x_4],[x__1(1),x__1(2),x__1(3),x__1(4)]);
s1 = -subs(P_gradient_sym, [x_1, x_2,x_3,x_4],[x__1(1),x__1(2),x__1(3),x__1(4)]);
lamda_star1=(transpose(g_f1)*g_f1)/(transpose(s1)*A1*s1);
x__2=x__1+lamda_star1*s1;e=1;
g_f2 = subs(P_gradient_sym, [x_1, x_2,x_3,x_4],[x__2(1),x__2(2),x__2(3),x__2(4)]);
grad_norm=norm(g_f2,inf);
k=2;
while (grad_norm>eP)
    A2 = subs(P_hessian_sym, [x_1, x_2,x_3,x_4],[x__2(1),x__2(2),x__2(3),x__2(4)]);
    beta2=(transpose(g_f2)*g_f2)/((transpose(g_f1)*g_f1));
    s2=-g_f2+beta2*s1;
    lamda_star2=(transpose(g_f2)*g_f2)/(transpose(s2)*A2*s2);
    x__3=x__2+lamda_star2*s2;
    x__2=x__3;
    g_f2 = subs(P_gradient_sym, [x_1, x_2,x_3,x_4],[x__2(1),x__2(2),x__2(3),x__2(4)]);
    grad_norm_prev=grad_norm;
    grad_norm=norm(g_f2, Inf);
    k=k+1;
    if (abs(grad_norm-grad_norm_prev)<e)
        break;
    end 
end 
fprintf('---------------The Powell’s quartic function using Reeves CG Method----------------\n');
fprintf('The final value at which function is min [%f,%f,%f,%f]\n',x__2(1),x__2(2),x__2(3),x__2(4));
end 