% Quadratic Interpolation Method
% Inputs:
%   t0: Initial guess for the minimum value
%   e1: Tolerance for the convergence criterion
%   f: Objective function to be minimized

function Quad(t0,e1,f)
A=0;
FA=f(A);F1=f(t0);
if F1>FA
    FC=F1;
    FB=f(t0/2);
    t=(t0/2);
else 
    while true
    FB=F1;
    F2=f(2*t0);
    if F1<F2
        FC=F2;
        t=t0;
        break;
    else 
        F1=F2;
        t0=2*t0;
        t=t0;
    end 
    end 
end 
A=0;B=t;C=2*t;
a=FA;
b=(4*FB-3*FA-FC)/(2*t);
c=(FC+FA-2*FB)/(2*(t^2));
h=@(lamda) a+b*(lamda)+c*(lamda^2);
lamda_star=(-b)/(2*c);
throd=abs((h(lamda_star)-f(lamda_star))/f(lamda_star));
k=1;
while (throd>=e1)
    if(A>lamda_star)
        FC=FB;C=B;
        FB=FA;B=A;
        FA=f(lamda_star);A=lamda_star;
    elseif (A<lamda_star && lamda_star>B )
        FC=FB;C=B;
        FB=f(lamda_star);B=lamda_star;
    elseif (lamda_star>B && C<lamda_star) 
        FA=FB;A=B;
        FB=f(lamda_star);B=lamda_star;
    else 
        FA=FB;A=B;
        FB=FC;B=C;
        FC=f(lamda_star);C=lamda_star;
    end
    break_throd=throd;
    a=FA;
    b=(4*FB-3*FA-FC)/(2*t);
    c=(FC+FA-2*FB)/(2*(t^2));
    h=@(lamda) a+b*(lamda)+c*(lamda^2);
    lamda_star=(-b)/(2*c);
    throd=abs((h(lamda_star)-f(lamda_star))/f(lamda_star));
    if (throd>=break_throd)
        break;
    end 
    k=k+1;
end 
fprintf('--------------------The Quadratic interpolation Method--------------------\n');
fprintf('the min value is %f and it functon value is %f\n',lamda_star,f(lamda_star));
end 