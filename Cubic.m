% Cubic Interpolation Method
% Inputs:
%   t0: Initial guess for the minimum value
%   e1: Tolerance for the convergence criterion
%   f: Objective function to be minimized

function Cubic(t0,e1,f)
df_dx = diff(f(sym('x')));
f_das = matlabFunction(df_dx);
A=0;
while f_das(A)>0
    if (f_das(A-t0)<f_das(A))
        A=A-t0;
    else
        A=A+t0;
    end 
end 
B=0;
while f_das(B)<=0
    if (f_das(B+t0)>f_das(B))
        B=B+t0;
    else
        B=B-t0;
    end 
end 
k=1;
while true
    FA=f(A);FA_das=f_das(A);
    FB=f(B);FB_das=f_das(B);
    Z = 3*(FA-FB)/(B-A)+ FA_das+ FB_das;
    AB_Diff_pow2=(A-B)^2;
    d = (2*Z + FA_das + FB_das)/(3*AB_Diff_pow2);
    c  = -((A+B)*Z + B*FA_das + A*FB_das)/AB_Diff_pow2;
    b  =  (B^2*FA_das+A^2*FB_das+2*A*B*Z)/AB_Diff_pow2;
    a = FA - b*A - c*A^2 - d*A^3;
    h=@(lamda) a+b*(lamda)+c*(lamda^2)+d*(lamda^3);
    lamda_star1=(-c+sqrt((c^2-3*b*d)))/(3*d);
    lamda_star2=(-c-sqrt((c^2-3*b*d)))/(3*d);
    if (A<=lamda_star1 && lamda_star1<=B)
        lamda_star=lamda_star1;
    else 
        lamda_star=lamda_star2;
    end 
    throd1=abs((h(lamda_star)-f(lamda_star))/f(lamda_star));
    throd2=abs(f_das(lamda_star));
    if (throd1<=e1 ||throd2<=e1)
        break;
    end 
        if (f_das(lamda_star)<0)
            A=lamda_star;
        else 
            B=lamda_star;
        end
  k=k+1;

end
fprintf('--------------------The Cubic interpolation Method--------------------\n');
fprintf('the min value is %f , it is functon value is %f and it is function dervative %f \n',lamda_star,f(lamda_star),f_das(lamda_star));
end 